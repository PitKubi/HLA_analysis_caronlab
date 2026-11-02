#!/usr/bin/env python3
"""Re-emit the collaborator's HLA-I prediction CSVs with the 9-mer binding core.

The collaborator ships one CSV per PTM in the HLA-I layout
    Peptide,Best_Allele,EL_Rank,Binder,Alleles,Tool
and the HLA-II layout of the same pipeline inserts Core after Best_Allele.
HLA-I never carried a core, so this script takes the binding block of every
PTM sheet in the results workbook, attaches the NetMHCpan 4.1b core predicted
for that peptide under the Best_Allele the workbook records, and writes the
HLA-II column order.

The cores are NOT recomputed here: they are read from the NetMHCpan output
already produced by hla1_netmhcpan_cores.py, which scored every peptide
against the three JY alleles.

Usage (from the repository root):
  python3 hla1_export_cores_for_collaborator.py
"""

import csv
from collections import Counter, OrderedDict
from pathlib import Path

WORKBOOK = Path("HLA-I_JY_DDA_PTMs_ResultsSummary_01262026_3.xlsx")
NETMHCPAN_OUTPUT_DIR = Path("/tmp/hla1_netmhcpan_cores")
EXPORT_DIR = Path("export_hla1_cores_2026-11-02")
PER_PTM_FILENAME = "sequentially_NetMHCpan_core.csv"
COMBINED_FILENAME = "all_PTMs_combined.csv"
CORE_DETAIL_FILENAME = "cores_full_detail.csv"

EXPORT_COLUMNS = ["Peptide", "Best_Allele", "Core", "EL_Rank", "Binder",
                  "Alleles", "Tool"]
CORE_DETAIL_COLUMNS = ["Peptide", "Best_Allele", "Core", "Of", "Gp", "Gl",
                       "Ip", "Il"]

CORE_LENGTH = 9
SHORTEST_PEPTIDE_WITH_RELIABLE_CORE = 9
LONGEST_PREDICTED_PEPTIDE = 14

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
NONSTANDARD_RESIDUE_PLACEHOLDER = "X"

# Column indices in a NetMHCpan 4.1 `-p` result line:
# Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL
ALLELE_FIELD = 1
PEPTIDE_FIELD = 2
CORE_FIELD = 3
OFFSET_FIELD = 4
GAP_POSITION_FIELD = 5
GAP_LENGTH_FIELD = 6
INSERTION_POSITION_FIELD = 7
INSERTION_LENGTH_FIELD = 8
IDENTITY_FIELD = 10
MINIMUM_RESULT_FIELDS = 13
PEPTIDE_LIST_IDENTITY = "PEPLIST"


class BindingBlock:
    """One binding table inside one workbook sheet, plus where it exports to.

    `folder` is the collaborator's own HLA-I folder name. `row_filter` keeps
    only part of the block: the 2026-01-26 workbook merged what used to be two
    acetylation blocks into one, so the K / N-term / other split now has to be
    read back out of the localisation column.
    """

    def __init__(self, folder, sheet, peptide_column, allele_column,
                 rank_column, binder_column, row_filter=None):
        self.folder = folder
        self.sheet = sheet
        self.peptide_column = peptide_column
        self.allele_column = allele_column
        self.rank_column = rank_column
        self.binder_column = binder_column
        self.row_filter = row_filter


ACETYL_LOCALISATION_COLUMN = "Localized acetytilation (residue)"


def acetylation_is_n_terminal(localisation):
    return "N-term" in localisation


def acetylation_carries_a_lysine_site(localisation):
    return any(token.strip().endswith("K") and token.strip()[:-1].isdigit()
               for token in localisation.split(","))


def acetylation_is_neither_lysine_nor_n_terminal(localisation):
    return not (acetylation_is_n_terminal(localisation)
                or acetylation_carries_a_lysine_site(localisation))


# The peptide/EL/allele/binder column names are the ones data_loader.R already
# verified against this workbook. Alleles and Tool are not named there; they
# sit in the two columns immediately after Binder in the same block, which
# read_binding_block checks before using.
BINDING_BLOCKS = [
    BindingBlock("Cyst_C", "Cyst.",
                 "Peptide", "Best_Allele...8", "EL_Rank", "Binder"),
    BindingBlock("Deamid_N", "Deamid. (NQ)",
                 "Peptide...7", "Best_Allele...8", "EL_Rank...9", "Binder...10"),
    BindingBlock("Deamid_Q", "Deamid. (NQ)",
                 "Peptide...16", "Best_Allele...17", "EL_Rank...18", "Binder...19"),
    BindingBlock("Phospho", "Phospho.",
                 "Peptide...23", "Best_Allele...24", "EL_Rank", "Binder"),
    BindingBlock("Acetyl_K", "Acetyl.",
                 "Peptide...43", "Best_Allele...47", "EL_Rank", "Binder",
                 row_filter=(ACETYL_LOCALISATION_COLUMN,
                             acetylation_carries_a_lysine_site)),
    BindingBlock("Acetyl_N-term", "Acetyl.",
                 "Peptide...43", "Best_Allele...47", "EL_Rank", "Binder",
                 row_filter=(ACETYL_LOCALISATION_COLUMN,
                             acetylation_is_n_terminal)),
    BindingBlock("Acetyl_others", "Acetyl.",
                 "Peptide...43", "Best_Allele...47", "EL_Rank", "Binder",
                 row_filter=(ACETYL_LOCALISATION_COLUMN,
                             acetylation_is_neither_lysine_nor_n_terminal)),
    BindingBlock("Methyl_KR", "Methyl.",
                 "Peptide...29", "Best_Allele...30", "EL_Rank...31", "Binder...32"),
    BindingBlock("Methyl_others", "Methyl.",
                 "Peptide...20", "Best_Allele...21", "EL_Rank...22", "Binder...23"),
    BindingBlock("Dimethyl", "Dimethyl.",
                 "Peptide...11", "Best_Allele...12", "EL_Rank", "Binder"),
    BindingBlock("Citrullination_R", "Citrullination",
                 "Peptide...11", "Best_Allele...12", "EL_Rank", "Binder"),
    BindingBlock("Carbamid_C", "Carbamid.",
                 "Peptide...11", "Best_Allele...12", "EL_Rank", "Binder"),
    BindingBlock("N-glyco", "N-glyco",
                 "StrippedSequences...18", "Best_Allele...20", "EL_Rank", "Binder"),
    BindingBlock("SUMO", "SUMO",
                 "Peptide...31", "Best_Allele...32", "EL_Rank", "Binder"),
    BindingBlock("Ubiq-G", "G-Ubiq.",
                 "Peptide...13", "Best_Allele...14", "EL_Rank", "Binder"),
    BindingBlock("Ubiq-GG", "GG-Ubiq.",
                 "Peptide", "Best_Allele...20", "EL_Rank", "Binder"),
    BindingBlock("Oxidation_bio", "bioOxid.",
                 "Peptide...142", "Best_Allele...145", "EL_Rank", "Binder"),
    BindingBlock("Oxid-MWHF_artefact", "artOxid.",
                 "Peptide...36", "Best_Allele...39", "EL_Rank", "Binder"),
    BindingBlock("Nitrosylation", "Nitrosylation",
                 "Peptide...22", "Best_Allele...23", "EL_Rank", "Binder"),
    BindingBlock("Background", "Background (wo PTMs)",
                 "Peptide", "Best_Allele...8", "EL_Rank", "Binder"),
]

EXPECTED_TOOL_NAMES = {"NetMHCpan", "MHCflurry"}
JY_SHORT_ALLELES = {"A0201", "B0702", "C0702"}


def repair_duplicate_names_like_readxl(header_cells):
    """Name the columns the way readxl's default name repair does.

    data_loader.R addresses these blocks by names such as "Peptide...43", so
    the export has to resolve the same strings to the same columns.
    """
    names = []
    for position, value in enumerate(header_cells, start=1):
        text = "" if value is None else str(value).strip()
        names.append(text if text else f"...{position}")
    occurrences = Counter(names)
    return [f"{name}...{position}" if occurrences[name] > 1 else name
            for position, name in enumerate(names, start=1)]


def load_sheet(workbook, sheet_name):
    worksheet = workbook[sheet_name]
    row_iterator = worksheet.iter_rows(values_only=True)
    column_names = repair_duplicate_names_like_readxl(next(row_iterator))
    return column_names, list(row_iterator)


def cell_text(row, index):
    if index is None or index >= len(row) or row[index] is None:
        return ""
    return str(row[index]).strip()


# The Acetyl. sheet stores the "no binding allele" case as a numeric 0 rather
# than an empty cell; every other sheet leaves it blank and the collaborator's
# own CSVs render it as an empty field. Carrying the 0 through would invent an
# allele named "0", so it is normalised back to blank and counted.
EMPTY_ALLELE_LIST_PLACEHOLDER = "0"


def normalise_empty_allele_list(alleles_text):
    return "" if alleles_text == EMPTY_ALLELE_LIST_PLACEHOLDER else alleles_text


def format_el_rank(value):
    """Render the numeric EL_Rank the way the collaborator's CSVs render it.

    repr() of the float is the shortest string that reads back as the same
    double, so the exported text round-trips to the workbook number exactly.
    """
    return repr(float(value))


def read_binding_block(block, column_names, rows):
    """Pull one block's peptide rows out of an already-loaded sheet."""
    peptide_index = column_names.index(block.peptide_column)
    allele_index = column_names.index(block.allele_column)
    rank_index = column_names.index(block.rank_column)
    binder_index = column_names.index(block.binder_column)
    alleles_index = binder_index + 1
    tool_index = binder_index + 2

    filter_index = None
    filter_test = None
    if block.row_filter is not None:
        filter_column, filter_test = block.row_filter
        filter_index = column_names.index(filter_column)

    extracted = []
    for row in rows:
        peptide = cell_text(row, peptide_index)
        if not peptide:
            continue
        if filter_test is not None and not filter_test(cell_text(row, filter_index)):
            continue
        extracted.append({
            "Peptide": peptide,
            "Best_Allele": cell_text(row, allele_index),
            "EL_Rank": format_el_rank(row[rank_index]),
            "Binder": cell_text(row, binder_index),
            "Alleles": normalise_empty_allele_list(cell_text(row, alleles_index)),
            "Tool": cell_text(row, tool_index),
        })
    return extracted


def assert_block_columns_are_the_expected_ones(block, extracted):
    """Guard the positional Alleles/Tool lookup and the named allele lookup."""
    unexpected_tools = {row["Tool"] for row in extracted} - EXPECTED_TOOL_NAMES
    assert not unexpected_tools, (
        f"{block.folder}: column after Binder+1 holds {unexpected_tools}, "
        f"which are not predictor names")
    allele_tokens = {token for row in extracted
                     for token in row["Alleles"].split(",") if token}
    assert allele_tokens <= JY_SHORT_ALLELES, (
        f"{block.folder}: Alleles column holds {allele_tokens - JY_SHORT_ALLELES}")
    best_alleles = {row["Best_Allele"] for row in extracted}
    assert best_alleles <= JY_SHORT_ALLELES, (
        f"{block.folder}: Best_Allele holds {best_alleles - JY_SHORT_ALLELES}")


def encode_nonstandard_residues(peptide):
    return "".join(residue if residue in STANDARD_AMINO_ACIDS
                   else NONSTANDARD_RESIDUE_PLACEHOLDER
                   for residue in peptide)


def short_allele_name(netmhcpan_allele):
    return netmhcpan_allele.replace("HLA-", "").replace("*", "").replace(":", "")


def read_core_predictions(output_dir):
    """Map (peptide, allele) to the core and its alignment offsets."""
    predictions = {}
    for result_path in sorted(output_dir.glob("*.txt")):
        with open(result_path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < MINIMUM_RESULT_FIELDS:
                    continue
                if fields[IDENTITY_FIELD] != PEPTIDE_LIST_IDENTITY:
                    continue
                key = (fields[PEPTIDE_FIELD],
                       short_allele_name(fields[ALLELE_FIELD]))
                predictions[key] = {
                    "Core": fields[CORE_FIELD],
                    "Of": int(fields[OFFSET_FIELD]),
                    "Gp": int(fields[GAP_POSITION_FIELD]),
                    "Gl": int(fields[GAP_LENGTH_FIELD]),
                    "Ip": int(fields[INSERTION_POSITION_FIELD]),
                    "Il": int(fields[INSERTION_LENGTH_FIELD]),
                }
    return predictions


def has_reliable_core(peptide):
    """8-mers are excluded: NetMHCpan pads them to nine with a '-' wildcard."""
    return (SHORTEST_PEPTIDE_WITH_RELIABLE_CORE <= len(peptide)
            <= LONGEST_PREDICTED_PEPTIDE)


def attach_cores(extracted, predictions):
    """Add the core predicted for each peptide under its Best_Allele."""
    exported = []
    dropped_eight_mers = 0
    for row in extracted:
        if not has_reliable_core(row["Peptide"]):
            dropped_eight_mers += 1
            continue
        prediction = predictions[(encode_nonstandard_residues(row["Peptide"]),
                                  row["Best_Allele"])]
        exported.append(dict(row, **prediction))
    return exported, dropped_eight_mers


def write_per_ptm_csv(folder, exported, export_dir):
    directory = export_dir / folder
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / PER_PTM_FILENAME
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPORT_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(exported)
    return path


def write_combined_csv(exported_by_folder, export_dir):
    path = export_dir / COMBINED_FILENAME
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["PTM"] + EXPORT_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        for folder, exported in exported_by_folder.items():
            for row in exported:
                writer.writerow(dict(row, PTM=folder))
    return path


def write_core_detail_csv(exported_by_folder, export_dir):
    """One row per distinct peptide/allele pair across the whole export."""
    detail = OrderedDict()
    for exported in exported_by_folder.values():
        for row in exported:
            detail.setdefault((row["Peptide"], row["Best_Allele"]), row)
    path = export_dir / CORE_DETAIL_FILENAME
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CORE_DETAIL_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in detail.values():
            writer.writerow(row)
    return path, len(detail)


def build_export():
    import openpyxl

    predictions = read_core_predictions(NETMHCPAN_OUTPUT_DIR)
    print(f"Loaded {len(predictions)} peptide/allele core predictions")

    workbook = openpyxl.load_workbook(WORKBOOK, read_only=True, data_only=True)
    loaded_sheets = {}
    exported_by_folder = OrderedDict()
    counts = []

    for block in BINDING_BLOCKS:
        if block.sheet not in loaded_sheets:
            loaded_sheets[block.sheet] = load_sheet(workbook, block.sheet)
        column_names, rows = loaded_sheets[block.sheet]
        extracted = read_binding_block(block, column_names, rows)
        assert_block_columns_are_the_expected_ones(block, extracted)
        exported, dropped = attach_cores(extracted, predictions)
        exported_by_folder[block.folder] = exported
        counts.append((block.folder, len(extracted), dropped, len(exported)))
        write_per_ptm_csv(block.folder, exported, EXPORT_DIR)

    write_combined_csv(exported_by_folder, EXPORT_DIR)
    detail_path, detail_rows = write_core_detail_csv(exported_by_folder, EXPORT_DIR)

    print(f"\n{'PTM':<22}{'workbook':>10}{'8-mers':>9}{'exported':>10}")
    for folder, workbook_rows, dropped, exported_rows in counts:
        print(f"{folder:<22}{workbook_rows:>10}{dropped:>9}{exported_rows:>10}")
    totals = [sum(column) for column in zip(*[row[1:] for row in counts])]
    print(f"{'TOTAL':<22}{totals[0]:>10}{totals[1]:>9}{totals[2]:>10}")
    print(f"\n{detail_path}: {detail_rows} distinct peptide/allele pairs")


if __name__ == "__main__":
    build_export()
