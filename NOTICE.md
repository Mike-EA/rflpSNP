# Third-party data notice

`rflpSNP` bundles a curated restriction enzyme dataset
(`data/restriction_enzymes.rda`) derived from:

- **REBASE**, the restriction enzyme database: Roberts RJ, Vincze T,
  Posfai J, Macelis D. "REBASE-a database for DNA restriction and
  modification: enzymes, genes and genomes." *Nucleic Acids Research*.
  Data used: emboss files, version 404 (2024). Free for academic use
  with attribution; see https://rebase.neb.com for terms.
- **Biopython**, `Bio.Restriction.Restriction_Dictionary`
  (https://github.com/biopython/biopython), used as the machine-readable
  intermediate source. Biopython License Agreement (BSD-3-style).

See `data-raw/PROVENANCE.md` for full curation methodology.

If `rflpSNP` is published (GitHub, CRAN, or an accompanying paper), cite
REBASE per the terms above.
