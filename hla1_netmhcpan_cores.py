#!/usr/bin/env python3
"""Build the HLA-I 9-mer binding-core tables that the core figures read.

The collaborator's HLA-I workbook ships an EL_Rank per peptide but no binding
core, while the HLA-II workbook ships cores. To give the manuscript an HLA-I
counterpart of the HLA-II core panels we re-run NetMHCpan 4.1b over every
unique HLA-I peptide for the three JY alleles and keep the predicted core.

Outputs, both written into figure_panels/:
  data_hla1_cores.csv       one row per peptide: core and its alignment offsets
  data_ptm_sites_core.csv   one row per modified site: its position inside the core

Usage (from the repository root):
  python3 hla1_netmhcpan_cores.py              run the predictor, then build
  python3 hla1_netmhcpan_cores.py --reuse      reuse an earlier run's output

--reuse reads the chunk output already present in --work-dir and skips the
predictor entirely; without it the work directory is rebuilt from scratch.
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter, OrderedDict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

NETMHCPAN_EXECUTABLE = Path("/home/pk/tools/netMHCpan-4.1/netMHCpan")
JY_ALLELES = "HLA-A02:01,HLA-B07:02,HLA-C07:02"
DEFAULT_WORK_DIR = Path("/tmp/hla1_netmhcpan_cores")
DEFAULT_CHUNK_COUNT = 16

FIGURE_PANEL_DIR = Path("figure_panels")
PTM_SITES_CSV = FIGURE_PANEL_DIR / "data_ptm_sites.csv"
BACKGROUND_CSV = FIGURE_PANEL_DIR / "data_background.csv"
BINDING_CSV = FIGURE_PANEL_DIR / "data_binding.csv"
CORE_TABLE_CSV = FIGURE_PANEL_DIR / "data_hla1_cores.csv"
PTM_SITE_CORE_CSV = FIGURE_PANEL_DIR / "data_ptm_sites_core.csv"

CORE_LENGTH = 9
SHORTEST_PREDICTED_PEPTIDE = 8
LONGEST_PREDICTED_PEPTIDE = 14

# NetMHCpan aligns an 8-mer to the 9-mer core by inserting a wildcard '-',
# so every 8-mer core carries a gap character at a position the predictor
# chose rather than a real residue. The collaborator asked for 8-mers to be
# left out of the core analysis for exactly that reason.
SHORTEST_PEPTIDE_WITH_RELIABLE_CORE = 9

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
NONSTANDARD_RESIDUE_PLACEHOLDER = "X"

# Column indices in a NetMHCpan 4.1 `-p` result line:
# Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL [BindLevel]
ALLELE_FIELD = 1
PEPTIDE_FIELD = 2
CORE_FIELD = 3
OFFSET_FIELD = 4
GAP_POSITION_FIELD = 5
GAP_LENGTH_FIELD = 6
INSERTION_POSITION_FIELD = 7
INSERTION_LENGTH_FIELD = 8
IDENTITY_FIELD = 10
RANK_FIELD = 12
MINIMUM_RESULT_FIELDS = 13
PEPTIDE_LIST_IDENTITY = "PEPLIST"


def reconstruct_core_from_alignment(peptide, offset, gap_position, gap_length,
                                    insertion_position, insertion_length):
    """Rebuild NetMHCpan's Core string from the peptide and its alignment."""
    if gap_length > 0:
        return (peptide[offset:offset + gap_position]
                + peptide[offset + gap_position + gap_length:
                          offset + CORE_LENGTH + gap_length])
    if insertion_length > 0:
        return (peptide[offset:offset + insertion_position]
                + "-" * insertion_length
                + peptide[offset + insertion_position:
                          offset + CORE_LENGTH - insertion_length])
    return peptide[offset:offset + CORE_LENGTH]


def core_position_of_peptide_site(site, offset, gap_position, gap_length,
                                  insertion_position, insertion_length):
    """Map a 1-based peptide position to a 1-based core position P1..P9.

    Returns None when the residue is not represented in the core: it either
    protrudes beyond the aligned window or sits in the segment NetMHCpan
    deleted to fit the peptide into nine positions.
    """
    offset_from_core_start = (site - 1) - offset
    if offset_from_core_start < 0:
        return None
    if gap_length > 0:
        if offset_from_core_start < gap_position:
            core_index = offset_from_core_start
        elif offset_from_core_start < gap_position + gap_length:
            return None
        else:
            core_index = offset_from_core_start - gap_length
    elif insertion_length > 0:
        core_index = (offset_from_core_start
                      if offset_from_core_start < insertion_position
                      else offset_from_core_start + insertion_length)
    else:
        core_index = offset_from_core_start
    if 0 <= core_index <= CORE_LENGTH - 1:
        return core_index + 1
    return None


def read_csv_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def is_predictable_length(peptide):
    return SHORTEST_PREDICTED_PEPTIDE <= len(peptide) <= LONGEST_PREDICTED_PEPTIDE


def encode_nonstandard_residues(peptide):
    """Replace residues NetMHCpan cannot score with its wildcard character."""
    return "".join(residue if residue in STANDARD_AMINO_ACIDS
                   else NONSTANDARD_RESIDUE_PLACEHOLDER
                   for residue in peptide)


def collect_peptides_to_predict(ptm_site_rows, background_rows):
    """Unique 8-14mers across the modified and background peptide tables."""
    peptides = {row["Peptide"] for row in ptm_site_rows}
    peptides.update(row["Peptide"] for row in background_rows)
    return sorted(peptide for peptide in peptides if is_predictable_length(peptide))


def build_encoded_peptide_index(peptides):
    """Map each NetMHCpan-safe sequence back to the sequence in the CSVs."""
    index = {}
    for peptide in peptides:
        encoded = encode_nonstandard_residues(peptide)
        if encoded in index and index[encoded] != peptide:
            raise ValueError(
                f"Wildcard encoding collapses {peptide!r} and {index[encoded]!r}")
        index[encoded] = peptide
    return index


def write_peptide_chunk_files(peptides, work_dir, chunk_count):
    work_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = []
    for chunk_index in range(chunk_count):
        chunk_peptides = peptides[chunk_index::chunk_count]
        if not chunk_peptides:
            continue
        chunk_path = work_dir / f"chunk{chunk_index:02d}.pep"
        chunk_path.write_text("\n".join(chunk_peptides) + "\n")
        chunk_paths.append(chunk_path)
    return chunk_paths


def predict_one_chunk(chunk_path):
    """Score one peptide chunk against the three JY alleles."""
    result_path = chunk_path.with_suffix(".txt")
    scratch_dir = chunk_path.parent / f"tmp_{chunk_path.stem}"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, TMPDIR=str(scratch_dir))
    with open(result_path, "w") as result_handle:
        subprocess.run(
            [str(NETMHCPAN_EXECUTABLE), "-p", str(chunk_path), "-a", JY_ALLELES],
            check=True, stdout=result_handle, env=environment,
        )
    return result_path


def predict_all_chunks_in_parallel(chunk_paths):
    with ThreadPoolExecutor(max_workers=len(chunk_paths)) as pool:
        return list(pool.map(predict_one_chunk, chunk_paths))


def locate_existing_results(work_dir):
    result_paths = sorted(work_dir.glob("*.txt"))
    if not result_paths:
        raise FileNotFoundError(
            f"--reuse was given but {work_dir} holds no NetMHCpan output")
    return result_paths


def short_allele_name(netmhcpan_allele):
    """HLA-A*02:01 -> A0201, matching the Allele column of the shipped CSVs."""
    return netmhcpan_allele.replace("HLA-", "").replace("*", "").replace(":", "")


def parse_prediction_records(result_paths):
    """Read every peptide/allele prediction out of the NetMHCpan output."""
    records = []
    for result_path in result_paths:
        with open(result_path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < MINIMUM_RESULT_FIELDS:
                    continue
                if fields[IDENTITY_FIELD] != PEPTIDE_LIST_IDENTITY:
                    continue
                records.append({
                    "encoded_peptide": fields[PEPTIDE_FIELD],
                    "allele": short_allele_name(fields[ALLELE_FIELD]),
                    "core": fields[CORE_FIELD],
                    "offset": int(fields[OFFSET_FIELD]),
                    "gap_position": int(fields[GAP_POSITION_FIELD]),
                    "gap_length": int(fields[GAP_LENGTH_FIELD]),
                    "insertion_position": int(fields[INSERTION_POSITION_FIELD]),
                    "insertion_length": int(fields[INSERTION_LENGTH_FIELD]),
                    "rank": float(fields[RANK_FIELD]),
                })
    return records


def assert_cores_reconstruct(records):
    """The alignment offsets must regenerate every Core the predictor printed.

    CoreSite is derived from those same offsets, so a single mismatch would
    silently shift modified residues to the wrong core position.
    """
    mismatches = [record for record in records
                  if reconstruct_core_from_alignment(
                      record["encoded_peptide"], record["offset"],
                      record["gap_position"], record["gap_length"],
                      record["insertion_position"],
                      record["insertion_length"]) != record["core"]]
    assert not mismatches, (
        f"core reconstruction failed for {len(mismatches)} of {len(records)} "
        f"predictions, first: {mismatches[0]}")
    print(f"Core reconstruction check: {len(records)}/{len(records)} exact")


def build_preferred_allele_index(binding_rows, background_rows):
    """The allele each peptide was already reported under by the collaborator.

    Anchoring to their Best_Allele keeps the core figures on the same
    peptide/allele pairing as every other HLA-I panel. The binding table wins
    over the background table, and its first occurrence wins within itself,
    because a peptide carrying several PTMs can be listed more than once.
    """
    preferred = {}
    for row in binding_rows:
        preferred.setdefault(row["Peptide"], row["Allele"])
    for row in background_rows:
        preferred.setdefault(row["Peptide"], row["Allele"])
    return preferred


def select_record_per_peptide(records, encoded_peptide_index, preferred_alleles):
    """Keep one prediction per peptide: its reported allele, else the best rank."""
    by_peptide = {}
    for record in records:
        peptide = encoded_peptide_index[record["encoded_peptide"]]
        by_peptide.setdefault(peptide, []).append(record)

    selected = OrderedDict()
    for peptide in sorted(by_peptide):
        candidates = by_peptide[peptide]
        preferred_allele = preferred_alleles.get(peptide)
        match = next((record for record in candidates
                      if record["allele"] == preferred_allele), None)
        if match is not None:
            selected[peptide] = (match, "Moh_Best_Allele")
        else:
            selected[peptide] = (min(candidates, key=lambda r: r["rank"]),
                                 "lowest_rank")
    return selected


def has_reliable_core(peptide):
    return (SHORTEST_PEPTIDE_WITH_RELIABLE_CORE <= len(peptide)
            <= LONGEST_PREDICTED_PEPTIDE)


def write_core_table(selected_records, path):
    rows_written = 0
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Peptide", "Allele", "Core", "Of", "Gp", "Gl", "Ip",
                         "Il", "EL_Rank", "AlleleSource"])
        for peptide, (record, allele_source) in selected_records.items():
            if not has_reliable_core(peptide):
                continue
            writer.writerow([peptide, record["allele"], record["core"],
                             record["offset"], record["gap_position"],
                             record["gap_length"], record["insertion_position"],
                             record["insertion_length"], record["rank"],
                             allele_source])
            rows_written += 1
    print(f"Wrote {path} ({rows_written} peptides)")


def build_ptm_site_core_rows(ptm_site_rows, selected_records):
    core_rows = []
    for row in ptm_site_rows:
        peptide = row["Peptide"]
        if not has_reliable_core(peptide) or peptide not in selected_records:
            continue
        record, _ = selected_records[peptide]
        site = row["Site"]
        core_site = None
        if site not in ("", "NA"):
            core_site = core_position_of_peptide_site(
                int(float(site)), record["offset"], record["gap_position"],
                record["gap_length"], record["insertion_position"],
                record["insertion_length"])
        core_rows.append({
            "PTM": row["PTM"],
            "Peptide": peptide,
            "Site": site,
            "Residue": row["Residue"],
            "Length": row["Length"],
            "Allele": record["allele"],
            "Core": record["core"],
            "CoreSite": "" if core_site is None else core_site,
        })
    return core_rows


def write_ptm_site_core_table(core_rows, path):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "PTM", "Peptide", "Site", "Residue", "Length", "Allele", "Core",
            "CoreSite"])
        writer.writeheader()
        writer.writerows(core_rows)
    print(f"Wrote {path} ({len(core_rows)} modified sites)")


def print_per_ptm_in_core_counts(ptm_site_rows, core_rows):
    """Report, per PTM, how much of the 9-14mer evidence lands inside a core."""
    eligible = Counter(row["PTM"] for row in ptm_site_rows
                       if has_reliable_core(row["Peptide"]))
    with_core = Counter(row["PTM"] for row in core_rows)
    in_core = Counter(row["PTM"] for row in core_rows if row["CoreSite"] != "")

    print("\nPer-PTM core coverage (9-14mer modified sites):")
    print(f"  {'PTM':<22} {'sites':>7} {'with core':>10} {'in P1-P9':>9} {'in-core %':>10}")
    for ptm, total in eligible.most_common():
        fraction = 100 * in_core[ptm] / total if total else 0.0
        print(f"  {ptm:<22} {total:>7} {with_core[ptm]:>10} {in_core[ptm]:>9} "
              f"{fraction:>9.1f}%")
    total_sites = sum(eligible.values())
    total_in_core = sum(in_core.values())
    print(f"  {'TOTAL':<22} {total_sites:>7} {sum(with_core.values()):>10} "
          f"{total_in_core:>9} {100 * total_in_core / total_sites:>9.1f}%")


def parse_command_line():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR,
                        help="directory holding the peptide chunks and the "
                             "NetMHCpan output")
    parser.add_argument("--reuse", action="store_true",
                        help="reuse the NetMHCpan output already in --work-dir")
    parser.add_argument("--chunks", type=int, default=DEFAULT_CHUNK_COUNT,
                        help="number of peptide chunks scored in parallel")
    return parser.parse_args()


def main():
    arguments = parse_command_line()

    ptm_site_rows = read_csv_rows(PTM_SITES_CSV)
    background_rows = read_csv_rows(BACKGROUND_CSV)
    binding_rows = read_csv_rows(BINDING_CSV)

    peptides = collect_peptides_to_predict(ptm_site_rows, background_rows)
    print(f"Unique {SHORTEST_PREDICTED_PEPTIDE}-{LONGEST_PREDICTED_PEPTIDE}mers "
          f"to score: {len(peptides)}")
    encoded_peptide_index = build_encoded_peptide_index(peptides)

    if arguments.reuse:
        result_paths = locate_existing_results(arguments.work_dir)
        print(f"Reusing {len(result_paths)} NetMHCpan result files from "
              f"{arguments.work_dir}")
    else:
        if not NETMHCPAN_EXECUTABLE.exists():
            sys.exit(f"NetMHCpan not found at {NETMHCPAN_EXECUTABLE}")
        chunk_paths = write_peptide_chunk_files(
            [encode_nonstandard_residues(peptide) for peptide in peptides],
            arguments.work_dir, arguments.chunks)
        print(f"Scoring {len(chunk_paths)} chunks against {JY_ALLELES}")
        result_paths = predict_all_chunks_in_parallel(chunk_paths)

    records = parse_prediction_records(result_paths)
    print(f"Parsed {len(records)} predictions")
    assert_cores_reconstruct(records)

    preferred_alleles = build_preferred_allele_index(binding_rows, background_rows)
    selected_records = select_record_per_peptide(
        records, encoded_peptide_index, preferred_alleles)

    write_core_table(selected_records, CORE_TABLE_CSV)
    core_rows = build_ptm_site_core_rows(ptm_site_rows, selected_records)
    write_ptm_site_core_table(core_rows, PTM_SITE_CORE_CSV)
    print_per_ptm_in_core_counts(ptm_site_rows, core_rows)


if __name__ == "__main__":
    main()
