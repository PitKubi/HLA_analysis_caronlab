HLA-I PTM peptides with NetMHCpan 9-mer binding cores
=====================================================
Prepared 2026-11-02 (Caron lab / Yale, PTM stability project).

These are the HLA-I per-PTM prediction tables with a Core column added, so
that the HLA-I data can be analysed at core positions P1..P9 the same way the
HLA-II data already is. The folder names, the file layout and the column order
follow the HLA-II export:

    Peptide,Best_Allele,Core,EL_Rank,Binder,Alleles,Tool

Core is the only column we added. Everything else is copied verbatim from the
binding block of the corresponding sheet of

    HLA-I_JY_DDA_PTMs_ResultsSummary_01262026_3.xlsx


Contents
--------
<PTM>/sequentially_NetMHCpan_core.csv
    One file per PTM, in the folder names used in HLA-I_PTMs_01062026/.
all_PTMs_combined.csv
    Every row of every per-PTM file, with a leading PTM column.
cores_full_detail.csv
    Peptide, Best_Allele, Core, Of, Gp, Gl, Ip, Il.
    One row per distinct peptide/allele pair. The Of/Gp/Gl/Ip/Il fields are
    needed to map a modification position into core coordinates (see below).


How the cores were produced
---------------------------
The HLA-I workbook reports EL_Rank but no core, so we ran the predictor
ourselves. NetMHCpan 4.1b (local install, /home/pk/tools/netMHCpan-4.1) was
run over every unique HLA-I peptide against the three JY class-I alleles:

    netMHCpan -p <peptide_list.pep> -a HLA-A02:01,HLA-B07:02,HLA-C07:02

with the peptide list split into 16 chunks scored in parallel. The Core, Of,
Gp, Gl, Ip and Il columns of the NetMHCpan `-p` output are taken as-is.

Each peptide is scored against all three alleles, but only one core is
exported per row: the core predicted under the allele that the workbook
already records in Best_Allele for that peptide in that sheet. That keeps the
core on exactly the same peptide/allele pairing as every other HLA-I number
in the manuscript. A peptide listed under different Best_Alleles in different
sheets therefore gets the core of its own sheet's allele.

Note that the EL_Rank/Binder values in the workbook come from NetMHCpan for
most rows and from MHCflurry for some (see the Tool column), whereas the Core
is always the NetMHCpan 4.1b core. MHCflurry does not report a core.


8-mers are excluded
-------------------
NetMHCpan aligns an 8-mer to the 9-position binding core by inserting a
wildcard '-' at a position the predictor picks, so an 8-mer "core" contains a
gap character rather than nine real residues and its P1..P9 coordinates are
not comparable with the rest. All 8-mers are therefore left out: 12,273 of the
90,755 workbook binding rows (13.5%). Only 9-14mers are exported, 78,482 rows.
Every exported Core is exactly 9 characters and contains no '-'.


Mapping a peptide position to a core position
---------------------------------------------
HLA-I cores are NOT always contiguous substrings of the peptide: for peptides
longer than 9 residues NetMHCpan deletes an internal stretch, and the core can
also start at an offset. Use Of/Gp/Gl/Ip/Il from cores_full_detail.csv.

Rebuilding the core from the peptide (all indices 0-based, Of = offset,
Gp/Gl = deletion position and length, Ip/Il = insertion position and length):

    if Gl > 0:   core = pep[Of : Of+Gp] + pep[Of+Gp+Gl : Of+9+Gl]
    elif Il > 0: core = pep[Of : Of+Ip] + '-'*Il + pep[Of+Ip : Of+9-Il]
    else:        core = pep[Of : Of+9]

Mapping a 1-based peptide position `site` to a 1-based core position P1..P9:

    d = (site - 1) - Of
    if d < 0:                       not in the core (before the aligned window)
    if Gl > 0:
        if d < Gp:                  core_index = d
        elif d < Gp + Gl:           not in the core (inside the deleted stretch)
        else:                       core_index = d - Gl
    elif Il > 0:
        core_index = d if d < Ip else d + Il
    else:
        core_index = d
    P = core_index + 1 if 0 <= core_index <= 8 else not in the core

A modified residue that falls outside the aligned window, or inside the
deleted stretch, has no core position and should be counted as "not in core"
rather than forced to a P number.

In this export Il is 0 everywhere (insertions only occur for 8-mers, which
are excluded), so in practice only the Gl and Of branches are ever taken:
37,744 of the 68,804 peptide/allele pairs have the core as the first nine
residues, 30,889 have an internal deletion (Gl > 0) and 781 start at an
offset (Of > 0). The Il branch is given above only for completeness.


Row counts per PTM
------------------
PTM folder              workbook rows   8-mers dropped   exported
Cyst_C                           2355              177       2178
Deamid_N                          958               14        944
Deamid_Q                          595               13        582
Phospho                           702               25        677
Acetyl_K                           57                8         49
Acetyl_N-term                     312               24        288
Acetyl_others                      23                0         23
Methyl_KR                          93                5         88
Methyl_others                     140               17        123
Dimethyl                          137                7        130
Citrullination_R                  138                0        138
Carbamid_C                         57                0         57
N-glyco                            55                2         53
SUMO                               30                2         28
Ubiq-G                            235                6        229
Ubiq-GG                           426               49        377
Oxidation_bio                     955               12        943
Oxid-MWHF_artefact               6382              380       6002
Nitrosylation                       5                1          4
Background                      77100            11531      65569
TOTAL                           90755            12273      78482


Three things that differ from the 2026-01-06 HLA-I folders
----------------------------------------------------------
1. Acetylation. The 2026-01-06 workbook had two separate acetylation binding
   blocks, which became the Acetyl_K and Acetyl_N-term folders. The
   2026-01-26 workbook merged them into one block of 381 rows carrying a
   "Localized acetytilation (residue)" column. We split that block back out:
   a row goes to Acetyl_N-term if its localisation string contains "N-term",
   and to Acetyl_K if it contains a lysine site such as "3K". Eleven rows
   carry both and appear in both files, exactly as they did in the 2026-01-06
   export. 23 rows localise acetylation only to S/T/C/H/Y and belong to
   neither of the two original folders; rather than drop them they are in an
   extra folder, Acetyl_others. Please tell us if you would prefer them
   merged elsewhere or dropped.

2. Nitrosylation is new in this workbook and had no folder before. Its binding
   block holds 5 rows (4 after the 8-mer filter). The folder is named
   Nitrosylation.

3. Repeated peptides. The bioOxid., artOxid., N-glyco and Acetyl. binding
   blocks list some peptides on more than one row. Those rows are kept, so the
   exported file has the same number of rows as the workbook block:
   Oxidation_bio 217 peptides on 328 extra rows, Oxid-MWHF_artefact 543
   peptides on 764 extra rows, N-glyco 5 peptides on 10 extra rows, Acetyl_K
   13 peptides on 13 extra rows. The 2026-01-06 export had collapsed these.

One data-cleaning step was applied to the values copied from the workbook: in
the Acetyl. sheet the "Alleles" cell of a non-binder holds a numeric 0 instead
of being left empty (319 rows). Every other sheet, and your own 2026-01-06
CSVs, leave it empty, so those cells are exported as an empty field.


Verification
------------
Checked against the workbook, re-read independently of any intermediate file:
  - all 78,482 Cores are 9 characters with no '-' wildcard;
  - all 68,804 distinct peptide/allele cores rebuild exactly from
    Of/Gp/Gl/Ip/Il with the rule above (100%);
  - Peptide, Best_Allele, Binder, Alleles and Tool match the workbook cell for
    cell in all 78,482 rows (392,410 field comparisons, 0 mismatches), and
    every exported EL_Rank reads back as the exact workbook number;
  - per-PTM row counts equal the workbook 9-14mer row counts for all 20
    folders;
  - no peptide repeats in a file unless it repeats in the workbook block.

For the 2026-01-06 folders, 78,430 of our rows have a peptide in common with
your file. All of them agree on EL_Rank. 47 differ in one field only: 41 are
Best_Allele tie-breaks on non-binders whose EL_Rank is a capped, tied value
(70, 75, 95, 100) and 6 are the same Alleles set written in a different order.

The export was produced by hla1_export_cores_for_collaborator.py in the
project repository; the cores come from hla1_netmhcpan_cores.py in the same
repository.
