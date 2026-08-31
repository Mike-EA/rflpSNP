# PCR-RFLP user guide

This guide describes a PCR-RFLP workflow in `rflpSNP`. It is intended for
known SNPs and a verifiable reference sequence. Read the [experimental design guide](GUIA_DISENO_EXPERIMENTAL.md)
and [limitations](LIMITACIONES_Y_VALIDACION.md) first.

## Requirements and inputs

Install `rflpSNP` and its dependencies as described in the [main README](../README.md).
Then check the Tm backend:

```r
library(rflpSNP)
packageVersion("TmCalculator")
check_tm_backend()
```

The validated version is `TmCalculator` 1.0.8. You need a FASTA sequence, a
flank immediately preceding the SNP under the `locate_snp()` convention, an
enzyme motif in IUPAC notation, and its cut position. HinfI, for example,
recognises `GANTC`; for `G^ANTC`, the cut offset is `1`.

```r
gene_seq <- read_gene_fasta("reference.fasta")
snp <- locate_snp(gene_seq, flank_seq = "5_PRIME_SEQUENCE_BEFORE_THE_SNP")
snp
```

Confirm that `snp$snp_base` is the expected reference allele. If the flank is
reported more than once, verify the coordinate with a longer flank or a more
specific reference before continuing.

## Step-by-step workflow

### 1. Generate and review candidates

```r
design <- design_primers(gene_seq, snp_pos = snp$snp_pos)
design$top
design$best_pair
```

`best_pair` is the leading candidate according to the internal score. Inspect
`top`, `n_candidates_forward`, `n_candidates_reverse`, and `n_pairs_total` as
well. The result does not establish enzyme discrimination or genome-wide primer
specificity. For difficult regions, change one parameter at a time and record
the reason; reasonable changes include the search window, amplicon range, Tm,
or GC limits.

### 2. Simulate the amplicon

```r
bp <- design$best_pair
pcr <- simulate_pcr(gene_seq, bp$forward_seq, bp$reverse_seq)
pcr
```

`n_sites_fwd` and `n_sites_rev` report potential matches in the FASTA used.
Values above one require independent specificity review.

### 3. Locate and review the restriction site

```r
site <- find_restriction_site(pcr, enzyme_motif = "GANTC", cut_offset = 1)
site$sites
```

Check both the number of sites and the cut position. A single site is not
sufficient by itself: it must differ between alleles and yield resolvable
fragments. Use `plot_amplicon_map()` and `plot_sequence_map()` to inspect the
layout.

### 4. Display expected bands and export the design

```r
simulate_gel(
  amplicon_size = pcr$size,
  fragment_sizes = c(146, 38),
  genotype_labels = c("ref/ref", "ref/alt", "alt/alt")
)
export_primers_txt(design, "primers_PCR_RFLP_results.txt")
```

Replace the example sizes with values derived from the confirmed cut. The plot
is an idealised representation, not a measurement of intensity or experimental
mobility. Archive the report with the FASTA, enzyme source, package version,
and Tm conditions.

## Quick pipeline

After reviewing each stage, run the workflow with one call:

```r
result <- run_pcr_rflp_pipeline(
  gene_seq,
  flank_seq = "5_PRIME_SEQUENCE_BEFORE_THE_SNP",
  enzyme_motif = "GANTC",
  cut_offset = 1,
  output_file = "primers_PCR_RFLP_results.txt"
)
```

The result contains `snp`, `design`, `pcr_result`, and `restriction_result`.
If the motif is absent from the reference allele, model and review the alternate
allele explicitly. See the [troubleshooting guide](CASOS_PROBLEMATICOS.md) and
complete the [validation checklist](LIMITACIONES_Y_VALIDACION.md) before
ordering oligonucleotides.
