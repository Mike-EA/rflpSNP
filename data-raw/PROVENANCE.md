# Provenance: `restriction_enzymes` dataset

## Source

The raw enzyme data was extracted from Biopython's `Bio.Restriction`
module, specifically the auto-generated file
`Bio/Restriction/Restriction_Dictionary.py`
(https://github.com/biopython/biopython), which is itself compiled from
**REBASE emboss files, version 404 (2024)** — the standard reference
database for restriction enzymes (Roberts, Vincze, Posfai & Macelis,
*Nucleic Acids Research*, REBASE database).

Biopython is distributed under the Biopython License Agreement / BSD-3
style license. REBASE data is freely redistributable for
academic/non-commercial use with attribution. Both sources are credited
here and should be credited in any publication or redistribution of
`rflpSNP` that relies on this dataset (see `NOTICE` file, to be added to
package root).

## Filtering criteria (catalog scope: "intermediate")

Starting universe: 1,088 enzymes in `rest_dict`.

Applied filters:

1. **`Commercially_available`** — enzyme is sold by at least one supplier
   (excludes purely computational/theoretical REBASE entries).
2. **`OneCut`** — enzyme has a single, well-defined cleavage position
   (excludes multi-cut / complex-mechanism enzymes not suited to simple
   RFLP fragment-size prediction).
3. **`site_length >= 4`** — excludes 2 edge-case entries (`FspEI`,
   `AbaSI`) whose REBASE "site" field is a 1-2 bp placeholder because
   their true specificity depends on cytosine methylation/
   hydroxymethylation state rather than DNA sequence. These do not fit
   a sequence-only in silico matching model and would produce massive
   false-positive hits if included.

**Result: 610 enzymes.**

### Important correction made during curation

REBASE/Biopython's `Defined` vs `Ambiguous` type tag was initially
assumed to mean "recognition site free of IUPAC ambiguity codes." This
assumption was tested and found **incorrect**: the tag reflects an
internal REBASE classification unrelated to sequence-level ambiguity or
to scientific certainty of the cut site. Using it as a filter would have
excluded widely-used, well-characterized enzymes such as `HinfI`
(already used in the package's original MTHFR C677T workflow), `AhdI`,
and `AarI`. The final filter (above) does not rely on this tag.

Recognition sites in the final dataset **do** legitimately contain IUPAC
ambiguity codes (e.g., `HinfI` = `GANTC`, `ApoI` = `RAATTY`) — this
reflects real enzyme biology. Sequence matching against these sites must
use `Biostrings::matchPattern(..., fixed = FALSE)`, consistent with how
`HinfI`/`GANTC` was already handled in the original scripts.

## Overhang sign convention (validated)

Biopython's `ovhg` field convention, confirmed against textbook
reference enzymes before use:

| Enzyme | Site   | ovhg | Overhang type | ovhg_seq |
|--------|--------|-----:|----------------|----------|
| EcoRI  | GAATTC |   -4 | 5' overhang    | AATT     |
| BamHI  | GGATCC |   -4 | 5' overhang    | GATC     |
| HindIII| AAGCTT |   -4 | 5' overhang    | AGCT     |
| PstI   | CTGCAG |   +4 | 3' overhang    | TGCA     |
| SmaI   | CCCGGG |    0 | blunt          | (none)   |
| HinfI  | GANTC  |   -3 | 5' overhang    | ANT      |

Convention used in `cut_type`: `ovhg > 0` → `3_overhang`,
`ovhg < 0` → `5_overhang`, `ovhg == 0` → `blunt`.

## Fields in `restriction_enzymes_source.csv`

| Column                  | Description |
|--------------------------|-------------|
| `enzyme`                 | Enzyme name (REBASE nomenclature) |
| `recognition_site`       | Recognition sequence, IUPAC codes allowed |
| `site_length`            | Length in bp of `recognition_site` |
| `fst5`                   | Cut position on top strand, relative to start of recognition site (REBASE convention) |
| `fst3`                   | Cut position on bottom strand, relative to end of recognition site |
| `ovhg`                   | Overhang length/sign (see convention above) |
| `ovhg_seq`               | Overhang sequence when applicable |
| `cut_type`               | `blunt` / `5_overhang` / `3_overhang` |
| `palindromic`            | Logical; FALSE for Type IIS-like non-palindromic recognition |
| `methylation_sensitive`  | Logical; enzyme activity is blocked/altered by genomic methylation (irrelevant for unmethylated PCR amplicons, kept as metadata only) |
| `n_suppliers`             | Count of commercial suppliers (REBASE) |
| `suppliers`               | Semicolon-separated supplier names |
| `rebase_id`               | REBASE database ID, for cross-reference |

## Regeneration

Run `data-raw/build_restriction_enzymes.R` to validate the CSV and
rebuild `data/restriction_enzymes.rda` and the roxygen documentation
stub. To regenerate the CSV itself from a newer REBASE/Biopython release,
re-run the extraction against an updated
`Restriction_Dictionary.py` and reapply the same three filters.
