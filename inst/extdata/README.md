# Enzyme database placeholder

Place your curated restriction enzyme CSV here as:

    restriction_enzymes.csv

with (at minimum) these columns, matching the schema used by
`read_enzyme_panel()`:

    enzyme,recognition_site,site_length,fst5,fst3,ovhg,ovhg_seq,cut_type,palindromic,methylation_sensitive,n_suppliers,suppliers,rebase_id

Once the file is in place, `pira_enzyme_panel()` (and therefore
`find_pira_sites()` / `design_pira_primers()`, which use it as their
default `enzyme_panel`) will load it automatically -- no code changes
needed.

To verify it loaded correctly after installing/reloading the package:

```r
p <- pira_enzyme_panel()
nrow(p)      # should match your CSV's enzyme count (after commercial_only filtering)
head(p)
```

If you see a warning about the bundled database not being found, the
file is either missing from this folder, named differently, or the
package needs to be reinstalled/reloaded (`devtools::load_all()` or
`devtools::install()`) after adding the file.

This file (`restriction_enzymes.csv`) should be deleted from this
directory before running `devtools::check()` if you have not yet
resolved its licensing status for redistribution -- see `../NOTICE.md`.
