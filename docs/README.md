# `rflpSNP` documentation

This directory contains the extended documentation for `rflpSNP`. The
[main README](../README.md) is intended for package installation and a first
analysis; these guides help users plan, review, and communicate a reproducible
genotyping assay.

## Reading paths

| If you need to… | Read |
|---|---|
| Understand the experimental principle and choose a strategy | [Experimental design guide](GUIA_DISENO_EXPERIMENTAL.md) |
| Run and review a PCR-RFLP analysis step by step | [PCR-RFLP user guide](PCR_RFLP_USER_GUIDE.md) |
| Design a tetra-primer assay for an SNP without a useful RFLP | [ARMS-PCR user guide](ARMS_PCR_USER_GUIDE.md) |
| Distinguish a computational prediction from experimental validation | [Limitations and validation](LIMITACIONES_Y_VALIDACION.md) |
| Investigate unexpected results or difficult inputs | [Troubleshooting guide](CASOS_PROBLEMATICOS.md) |
| Tune filters and review reproducible MTHFR/FTO examples | [Parameter tuning and worked examples](PARAMETER_TUNING_AND_WORKED_EXAMPLES.md) |
| Review the coordinate and orientation contract for ARMS-PCR | [ARMS-PCR conventions](ARMS_PCR_CONVENTIONS.md) |

The [ARMS-PCR development plan](ARMS_PCR_PLAN.md) records the historical
decisions and scope of that feature; it is not a user guide.

## Scope

`rflpSNP` generates and compares *in silico* candidates for known SNPs. It is
intended for teaching, assay planning, and reproducible analysis. It does not
confirm genomic specificity, PCR performance, enzyme digestion, or the
genotype of a biological sample. Every design must be reviewed and validated
under the applicable experimental and regulatory conditions.

## Shared conventions

- Primer sequences are shown 5'→3'.
- FASTA coordinates are 1-based and inclusive.
- The reference allele is the base in the FASTA sequence at the analysed
  coordinate; the alternate allele must be supplied in the same orientation.
- Amplicon and band sizes are expressed in base pairs (bp).
