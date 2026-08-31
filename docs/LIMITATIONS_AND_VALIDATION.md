# Limitations and Validation

`rflpSNP` produces *in silico* designs and predictions. A valid result is a
traceable assay hypothesis, not evidence that a reaction will work or that a
sample has a particular genotype.

## What the package covers

- SNP localisation from a flank in the supplied FASTA.
- Generation, filtering, and ranking of PCR-RFLP and tetra-primer ARMS-PCR
  primers.
- Tm calculation under declared conditions and heuristic screening of %GC,
  dimers, and hairpins.
- Product construction from the reference sequence and IUPAC motif searches.
- Amplicon diagrams and virtual band patterns.

## What the package does not cover

- Exhaustive specificity searches against genomes, transcriptomes, pseudogenes,
  paralogues, or population variants.
- Structural variation, CNVs, indels, or primer-binding-site polymorphisms,
  unless explicitly represented in a reviewed sequence.
- PCR efficiency, allelic bias, inhibition, contamination, DNA concentration,
  or other properties of a biological sample or thermocycler.
- Restriction-enzyme activity, purity, methylation sensitivity, buffer, time,
  or partial digestion.
- Actual gel resolution, band intensity, staining, loading, or small-fragment
  behaviour.
- Clinical validation, diagnosis, medical interpretation, or regulatory
  compliance.

## Important assumptions

`locate_snp()` assumes a relevant flank match exists; when several matches
occur, it uses the first and emits a warning. `simulate_pcr()` uses the first
appropriately oriented exact matches to build an amplicon. Its mismatch counts
are a local alert, not a genomic specificity analysis. `find_restriction_site()`
interprets the motif and `cut_offset` supplied by the user.

Tm values depend on salt, Mg2+, dNTP, and oligo concentrations as well as the
thermodynamic model. Dimer and hairpin estimates are based on complementarity
and prioritise candidates; they do not replace specialised thermodynamic
evaluation.

For ARMS-PCR, deliberate mismatches and 3' discrimination are computational
design rules. Selectivity also depends on polymerase, temperature,
concentrations, and template. A polymerase without 3'→5' proofreading activity
and experimental optimisation are required.

## Validation checklist before synthesis

1. Confirm the SNP, alleles, orientation, assembly, and reference version.
2. Check every primer against the relevant genomic or transcriptomic reference
   with an independent specificity tool.
3. Review Tm, self/heterodimers, hairpins, and 3' ends under actual PCR mixture
   conditions.
4. Confirm motif, cuts, additional sites, and digestion conditions in current
   enzyme documentation.
5. Ensure products can be distinguished with the selected agarose percentage,
   ladder, voltage, and detection system.
6. Define no-template, digestion, and known-genotype controls where possible;
   predefine how ambiguous patterns will be resolved.
7. Optimise annealing temperature, Mg2+, cycles, and—in ARMS-PCR—primer ratios
   and the second mismatch experimentally.
8. Confirm a subset of genotypes by an independent method before routine use.

## Responsible use

The package is intended for research, planning, and teaching. It must not be
the sole basis for clinical, diagnostic, or regulatory decisions. The assay
lead is responsible for analytical suitability and local requirements.
