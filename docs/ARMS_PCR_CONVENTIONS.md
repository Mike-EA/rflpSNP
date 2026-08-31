# ARMS-PCR coordinate and orientation conventions

This document defines the representation contract used by ARMS-PCR functions.

## Sequences and coordinates

- `gene_seq` is a `Biostrings::DNAString` reference sequence represented 5'→3'.
- Genomic coordinates (`snp_pos`, primer start and end) are 1-based, inclusive,
  and always refer to `gene_seq`.
- An interval from `start` to `end` has length `end - start + 1` and `start <= end`.
- `locate_snp()` returns `snp_pos` in this system; its `strand` field does not
  change the returned coordinate frame.

## Alleles and primers

- `ref_allele` is the uppercase base in `gene_seq[snp_pos]` and is validated
  before assay design.
- `alt_allele` is explicit, one of `A`, `C`, `G`, or `T`, uppercase, and differs
  from `ref_allele`.
- All primer sequences are stored and exported 5'→3'.
- A forward primer matches `gene_seq[start:end]`; a reverse primer is
  `reverseComplement(gene_seq[start:end])`.
- An allele-specific forward inner primer ends at the SNP; an allele-specific
  reverse inner primer begins at the SNP.
- Deliberate mismatches are described from the 3' end: `3'-2` and `3'-3` map to
  R indices `length(primer)-2` and `length(primer)-3`. The terminal SNP base is
  named `3'` and is not a deliberate mismatch.
