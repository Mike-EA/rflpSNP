# ARMS-PCR development plan

## Purpose

ARMS-PCR adds *in silico* tetra-primer assay design for biallelic SNPs when a
PCR-RFLP restriction-site change is not useful. The workflow produces four
primers, a text report, expected products for three genotypes, and a virtual
gel. It does not provide genomic specificity prediction, experimental
validation, PIRA-PCR, multiplexing, or multiallelic assays.

## Technical decision

Each assay uses two outer primers for a control product and two allele-specific
inner primers. Inner primers end at the SNP and carry a configurable deliberate
mismatch at 3'-2 or 3'-3. Candidate sets are preferred when diagnostic products
and the control product are distinguishable and all four primers have compatible
physicochemical properties with low heuristic dimer and hairpin risk.

## Implemented workflow

```text
Reference FASTA → SNP flank → reference coordinate and allele
                               + explicit alternate allele
                 → ARMS candidate design and ranking
                 → report → product simulation → virtual gel
```

The public interface consists of `design_arms_primers()`,
`export_arms_primers_txt()`, `simulate_arms_pcr()`, `simulate_arms_gel()`, and
`run_arms_pcr_pipeline()`.

## Quality expectations

The implementation is covered by synthetic, regression, orientation,
transition/transversion, filtering, and genotype-pattern tests. Documentation
and validation guidance emphasise that candidates require independent genomic
specificity assessment, specialised thermodynamic review, a non-proofreading
polymerase, and laboratory optimisation.

See the [user guide](ARMS_PCR_USER_GUIDE.md) for use and the
[conventions](ARMS_PCR_CONVENTIONS.md) for the coordinate contract.
