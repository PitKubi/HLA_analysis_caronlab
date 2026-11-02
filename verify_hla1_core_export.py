#!/usr/bin/env python3
"""Independent verification of export_hla1_cores_2026-11-02/.

Re-reads the source workbook from scratch (never the intermediate CSVs
produced for the figures) and checks every exported row against it.
"""

import csv
import random
import sys
from collections import Counter, OrderedDict
from pathlib import Path

REPO = Path("/home/pk/PycharmProjects/CaronLab_Yale/PTM_stability_caron")
sys.path.insert(0, str(REPO))

import openpyxl  # noqa: E402

from hla1_export_cores_for_collaborator import (  # noqa: E402
    BINDING_BLOCKS, EXPORT_DIR, PER_PTM_FILENAME, COMBINED_FILENAME,
    CORE_DETAIL_FILENAME, EXPORT_COLUMNS, CORE_DETAIL_COLUMNS, WORKBOOK,
    repair_duplicate_names_like_readxl, cell_text,
)

CORE_LENGTH = 9
failures = []


def check(label, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {label}{(' - ' + detail) if detail else ''}")
    if not ok:
        failures.append(label)


def read_csv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def rebuild_core(peptide, of, gp, gl, ip, il):
    if gl > 0:
        return peptide[of:of + gp] + peptide[of + gp + gl:of + CORE_LENGTH + gl]
    if il > 0:
        return (peptide[of:of + ip] + "-" * il
                + peptide[of + ip:of + CORE_LENGTH - il])
    return peptide[of:of + CORE_LENGTH]


# ---------------------------------------------------------------- load export
export_dir = REPO / EXPORT_DIR
per_ptm = OrderedDict()
for block in BINDING_BLOCKS:
    header, rows = read_csv(export_dir / block.folder / PER_PTM_FILENAME)
    assert header == EXPORT_COLUMNS, (block.folder, header)
    per_ptm[block.folder] = rows
check("header of every per-PTM CSV is exactly " + ",".join(EXPORT_COLUMNS), True,
      f"{len(per_ptm)} files")

combined_header, combined = read_csv(export_dir / COMBINED_FILENAME)
check("all_PTMs_combined.csv header", combined_header == ["PTM"] + EXPORT_COLUMNS)
detail_header, detail = read_csv(export_dir / CORE_DETAIL_FILENAME)
check("cores_full_detail.csv header", detail_header == CORE_DETAIL_COLUMNS)

all_rows = [row for rows in per_ptm.values() for row in rows]
check("combined row count == sum of per-PTM rows",
      len(combined) == len(all_rows), f"{len(combined)} vs {len(all_rows)}")

# ------------------------------------------------------------------- check (a)
bad_length = [r for r in all_rows if len(r["Core"]) != CORE_LENGTH]
bad_wildcard = [r for r in all_rows if "-" in r["Core"]]
check("(a) every Core is 9 characters",
      not bad_length, f"{len(all_rows)} rows checked, {len(bad_length)} bad")
check("(a) no Core contains a '-' wildcard",
      not bad_wildcard, f"{len(bad_wildcard)} bad")

# ------------------------------------------------------------------- check (b)
detail_by_key = {(r["Peptide"], r["Best_Allele"]): r for r in detail}
mismatch = []
for row in detail:
    rebuilt = rebuild_core(row["Peptide"], int(row["Of"]), int(row["Gp"]),
                           int(row["Gl"]), int(row["Ip"]), int(row["Il"]))
    if rebuilt != row["Core"]:
        mismatch.append((row, rebuilt))
check("(b) core reconstructs from Of/Gp/Gl/Ip/Il in cores_full_detail.csv",
      not mismatch, f"{len(detail)}/{len(detail)} exact"
      if not mismatch else f"{len(mismatch)} mismatches, first {mismatch[0]}")

missing_detail = [r for r in all_rows
                  if (r["Peptide"], r["Best_Allele"]) not in detail_by_key]
core_disagree = [r for r in all_rows
                 if (r["Peptide"], r["Best_Allele"]) in detail_by_key
                 and detail_by_key[(r["Peptide"], r["Best_Allele"])]["Core"] != r["Core"]]
check("(b) every exported peptide/allele pair is in cores_full_detail.csv",
      not missing_detail, f"{len(missing_detail)} missing")
check("(b) Core agrees between per-PTM files and cores_full_detail.csv",
      not core_disagree, f"{len(core_disagree)} disagreements")

exported_pairs = {(r["Peptide"], r["Best_Allele"]) for r in all_rows}
check("(b) cores_full_detail.csv has no rows beyond the export",
      len(detail) == len(exported_pairs), f"{len(detail)} vs {len(exported_pairs)}")

# ------------------------------------------------- checks (c) (d) (e) via xlsx
workbook = openpyxl.load_workbook(REPO / WORKBOOK, read_only=True, data_only=True)
loaded = {}
count_table = []
field_mismatches = []
rank_mismatches = []
duplicate_problems = []
spot_pool = []

for block in BINDING_BLOCKS:
    if block.sheet not in loaded:
        worksheet = workbook[block.sheet]
        iterator = worksheet.iter_rows(values_only=True)
        names = repair_duplicate_names_like_readxl(next(iterator))
        loaded[block.sheet] = (names, list(iterator))
    names, rows = loaded[block.sheet]
    ipep = names.index(block.peptide_column)
    iall = names.index(block.allele_column)
    irank = names.index(block.rank_column)
    ibind = names.index(block.binder_column)
    ifilter, ftest = (None, None)
    if block.row_filter is not None:
        ifilter = names.index(block.row_filter[0])
        ftest = block.row_filter[1]

    workbook_rows = []
    for excel_row_number, row in enumerate(rows, start=2):
        peptide = cell_text(row, ipep)
        if not peptide:
            continue
        if ftest is not None and not ftest(cell_text(row, ifilter)):
            continue
        workbook_rows.append((excel_row_number, peptide, row))

    nine_to_fourteen = [t for t in workbook_rows if 9 <= len(t[1]) <= 14]
    eight_mers = len(workbook_rows) - len(nine_to_fourteen)
    exported = per_ptm[block.folder]
    count_table.append((block.folder, block.sheet, block.peptide_column,
                        len(workbook_rows), eight_mers, len(nine_to_fourteen),
                        len(exported)))

    if len(exported) != len(nine_to_fourteen):
        continue

    for (excel_row_number, peptide, row), out in zip(nine_to_fourteen, exported):
        for column, index in (("Peptide", ipep), ("Best_Allele", iall),
                              ("Binder", ibind), ("Tool", ibind + 2)):
            if out[column] != cell_text(row, index):
                field_mismatches.append((block.folder, excel_row_number, column,
                                         out[column], cell_text(row, index)))
        workbook_alleles = cell_text(row, ibind + 1)
        if workbook_alleles == "0":
            workbook_alleles = ""
        if out["Alleles"] != workbook_alleles:
            field_mismatches.append((block.folder, excel_row_number, "Alleles",
                                     out["Alleles"], workbook_alleles))
        if float(out["EL_Rank"]) != float(row[irank]):
            rank_mismatches.append((block.folder, excel_row_number,
                                    out["EL_Rank"], row[irank]))
        spot_pool.append((block.folder, block.sheet, excel_row_number, out, row,
                          (ipep, iall, irank, ibind)))

    exported_counts = Counter(r["Peptide"] for r in exported)
    workbook_counts = Counter(p for _, p, _ in nine_to_fourteen)
    for peptide, n in exported_counts.items():
        if n != workbook_counts.get(peptide, 0):
            duplicate_problems.append((block.folder, peptide, n,
                                       workbook_counts.get(peptide, 0)))

check("(d) exported rows == workbook 9-14mer rows for every PTM",
      all(r[5] == r[6] for r in count_table),
      f"{sum(r[6] for r in count_table)} rows over {len(count_table)} PTMs")
check("(c) Peptide/Best_Allele/Binder/Alleles/Tool identical to the workbook",
      not field_mismatches,
      f"{len(spot_pool) * 5} field comparisons, {len(field_mismatches)} mismatches")
check("(c) EL_Rank round-trips to the workbook number exactly",
      not rank_mismatches,
      f"{len(spot_pool)} compared, {len(rank_mismatches)} mismatches")
check("(e) per-peptide multiplicity matches the workbook block",
      not duplicate_problems,
      f"{len(duplicate_problems)} files/peptides off"
      if duplicate_problems else "duplicates only where the workbook duplicates")

# ----------------------------------------------------------------- count table
print("\nPer-PTM row counts")
print(f"{'PTM folder':<22}{'sheet':<22}{'peptide column':<26}"
      f"{'workbook':>9}{'8-mers':>8}{'expected':>9}{'exported':>9}")
for folder, sheet, column, wb_rows, eights, expected, exported_n in count_table:
    print(f"{folder:<22}{sheet:<22}{column:<26}{wb_rows:>9}{eights:>8}"
          f"{expected:>9}{exported_n:>9}")
totals = [sum(c[i] for c in count_table) for i in (3, 4, 5, 6)]
print(f"{'TOTAL':<22}{'':<22}{'':<26}{totals[0]:>9}{totals[1]:>8}"
      f"{totals[2]:>9}{totals[3]:>9}")

# ------------------------------------------------------------- duplicate detail
dupes = [(folder, Counter(r["Peptide"] for r in rows))
         for folder, rows in per_ptm.items()]
print("\nFiles containing repeated peptides (workbook repeats them too)")
for folder, counter in dupes:
    repeated = {p: n for p, n in counter.items() if n > 1}
    if repeated:
        print(f"  {folder:<22} {len(repeated)} peptides repeated, "
              f"{sum(repeated.values()) - len(repeated)} extra rows")

# ------------------------------------------------------------------- check (f)
print("\n(f) 10 random rows re-read cell by cell from the workbook")
random.seed(20261102)
by_folder = OrderedDict()
for entry in spot_pool:
    by_folder.setdefault(entry[0], []).append(entry)
spot_sample = [random.choice(by_folder[folder])
               for folder in random.sample(list(by_folder), 10)]
print(f"{'PTM':<20}{'xlsx row':>9}  {'Peptide':<15}{'Allele':<8}{'Core':<11}"
      f"{'EL_Rank':<14}{'Binder':<11}{'Alleles':<18}{'Tool':<11}{'ok'}")
for folder, sheet, excel_row_number, out, row, (ipep, iall, irank, ibind) in \
        spot_sample:
    workbook_alleles = cell_text(row, ibind + 1)
    if workbook_alleles == "0":
        workbook_alleles = ""
    same = (out["Peptide"] == cell_text(row, ipep)
            and out["Best_Allele"] == cell_text(row, iall)
            and float(out["EL_Rank"]) == float(row[irank])
            and out["Binder"] == cell_text(row, ibind)
            and out["Alleles"] == workbook_alleles
            and out["Tool"] == cell_text(row, ibind + 2)
            and len(out["Core"]) == 9 and "-" not in out["Core"])
    detail_row = detail_by_key[(out["Peptide"], out["Best_Allele"])]
    rebuilt = rebuild_core(out["Peptide"], int(detail_row["Of"]),
                           int(detail_row["Gp"]), int(detail_row["Gl"]),
                           int(detail_row["Ip"]), int(detail_row["Il"]))
    same = same and rebuilt == out["Core"]
    print(f"{folder:<20}{excel_row_number:>9}  {out['Peptide']:<15}"
          f"{out['Best_Allele']:<8}{out['Core']:<11}{out['EL_Rank']:<14}"
          f"{out['Binder']:<11}{(out['Alleles'] or '(empty)'):<18}"
          f"{out['Tool']:<11}{'OK' if same else 'MISMATCH'}")
    print(f"{'  xlsx cells':<20}{'':>9}  {cell_text(row, ipep):<15}"
          f"{cell_text(row, iall):<8}{'':<11}{row[irank]!r:<14}"
          f"{cell_text(row, ibind):<11}{(workbook_alleles or '(empty)'):<18}"
          f"{cell_text(row, ibind + 2):<11}")
    if not same:
        failures.append(f"spot check {folder} row {excel_row_number}")

print("\n" + ("ALL CHECKS PASSED" if not failures
             else f"FAILURES: {failures}"))
