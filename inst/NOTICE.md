# Third-party data notice

## Restriction enzyme database (`inst/extdata/restriction_enzymes.csv`)

If this file is derived from REBASE (The Restriction Enzyme Database,
http://rebase.neb.com), note that REBASE data is distributed under a
**Creative Commons Attribution-NonCommercial (CC BY-NC)** license.

This means:

- Attribution to REBASE is required when redistributing the data
  (including as part of this R package).
- **Non-commercial use only.** If `rflpSNP` (or this data file
  specifically) is intended for any commercial use, this data cannot be
  redistributed under its current license without separate permission
  from REBASE.
- If `rflpSNP` is submitted to CRAN, Bioconductor, or JOSS, review
  whether a CC BY-NC dependency is compatible with that venue's policies
  before submission -- CRAN in particular is often used in commercial
  contexts, which can be a friction point for CC BY-NC-licensed bundled
  data. Consider either: (a) keeping the full CSV out of the package and
  distributing it separately (e.g. via `read_enzyme_panel()` pointing to
  a user-downloaded file, as before), or (b) confirming with REBASE that
  bundling is acceptable for your intended distribution channel.

Suggested attribution line (adjust to REBASE's current requested
citation, see http://rebase.neb.com/rebase/rebase.html):

> Restriction enzyme recognition sites and cut positions are derived from
> REBASE, The Restriction Enzyme Database (Roberts, R.J., Vincze, T.,
> Posfai, J., Macelis, D. (2023) Nucleic Acids Res. 43: D298-D299),
> http://rebase.neb.com. Used under CC BY-NC.

If this CSV is **not** REBASE-derived (e.g. hand-compiled from supplier
catalogs), delete this notice or replace it with the correct source and
license for your actual data.
