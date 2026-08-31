# Tetra-primer ARMS-PCR user guide

This guide describes the `rflpSNP` ARMS-PCR workflow for biallelic SNPs without
a useful restriction-site change. It returns a four-primer set, expected
products for three genotypes, a text report, a virtual gel, and nucleotide maps.

Read the [experimental design guide](GUIA_DISENO_EXPERIMENTAL.md) before using
this workflow. Install the package as described in the [main README](../README.md).

## Scope and limits

The module designs *in silico* candidates. It does not establish genomic
specificity, replace specialised thermodynamic validation, or guarantee
experimental performance. Use a polymerase without 3'→5' proofreading activity
and validate annealing temperature and deliberate mismatch in the laboratory.
Complete the [validation checklist](LIMITACIONES_Y_VALIDACION.md) before
synthesis.

## Basic workflow

```r
library(rflpSNP)

gene_seq <- read_gene_fasta("reference.fasta")
snp <- locate_snp(gene_seq, flank_seq = "5_PRIME_SEQUENCE_BEFORE_THE_SNP")

design <- design_arms_primers(
  gene_seq,
  snp_pos = snp$snp_pos,
  ref_allele = snp$snp_base,
  alt_allele = "A"
)

design$best_set
```

`ref_allele` must match the reference base at `snp_pos`. `alt_allele` must be
the other base in the same FASTA orientation.

## Interpreting a design

A valid result includes two outer primers for the control band, one reference-
specific inner primer, one alternate-specific inner primer, and three product
sizes: control, reference, and alternate. `design$diagnostics` summarises
individual filtering. If no valid set is found, inspect it together with
`design$exclusion_diagnostics`.

```r
design$diagnostics
design$exclusion_diagnostics
```

`max_raw_candidates_per_pool` limits candidates evaluated before the
thermodynamic calculation. It supports rapid exploration but does not replace
an exhaustive final search.

## Products, gel, and maps

```r
pcr <- simulate_arms_pcr(
  gene_seq, design$best_set,
  snp_pos = snp$snp_pos,
  ref_allele = snp$snp_base,
  alt_allele = "A"
)

pcr$products
simulate_arms_gel(pcr)
export_arms_primers_txt(design, "arms_design.txt")
export_arms_amplicon_map_txt(pcr, "arms_nucleotide_map.txt")
```

`ref/ref` shows control plus reference band; `ref/alt` shows all three bands;
and `alt/alt` shows control plus alternate band. The map records both strands,
FASTA coordinates, the SNP, and primers stored 5'→3'.

## Pipeline

```r
result <- run_arms_pcr_pipeline(
  gene_seq,
  flank_seq = "5_PRIME_SEQUENCE_BEFORE_THE_SNP",
  alt_allele = "A",
  output_file = "arms_design.txt"
)
```

The result retains `snp`, `design`, `report`, `pcr_result`, and `gel` for
reproducible review. Use the [troubleshooting guide](CASOS_PROBLEMATICOS.md)
when no set is found or the expected bands are ambiguous.
