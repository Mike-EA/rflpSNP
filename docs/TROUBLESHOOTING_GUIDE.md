# Troubleshooting Guide

Keep the error message, arguments, and FASTA provenance before changing a
parameter. Unrecorded changes make design attempts difficult to compare.

## The flank is not found or occurs more than once

**Possible causes:** wrong assembly, opposite orientation, a short flank,
incomplete sequence, or mixing transcript and genomic references.

**Actions:** confirm the reference version; use a longer, more specific flank;
check both orientations; and ensure the SNP lies within the FASTA. When several
matches occur, verify the coordinate independently rather than accepting the
first one automatically.

## The reference allele does not match the expected SNP

**Possible causes:** an assembly change, opposite-strand nomenclature, a
coordinate error, or a reference carrying an unexpected variant.

**Actions:** stop the design. Reconcile alleles and orientation before changing
`snp_offset` or primer criteria.

## No primers or valid pairs are generated

**Possible causes:** a short region, extreme GC content, repeats, restrictive
Tm or GC limits, incompatible amplicon sizes, or non-resolvable fragments.

**Actions:** identify the limiting criterion; extend the region when sequence
is available; change one criterion at a time and record it. If alternatives
remain poor, consider ARMS-PCR or another genotyping method.

## In-silico PCR reports several potential sites

**Interpretation:** additional matches were detected in the input FASTA at the
allowed mismatch level. This neither demonstrates nor excludes real
non-specific amplification.

**Actions:** perform a specificity search against the relevant genome, review
3' ends closely, and compare other candidates in `design$top`.

## The expected restriction site is absent

**Possible causes:** an incorrect motif or `cut_offset`, a reference allele
without the site, an orientation error, or an incorrect input sequence.

**Actions:** confirm the enzyme documentation and model the alternate allele.
Also check for multiple sites before interpreting a simple two-fragment pattern.

## Fragments are too close, too small, or not visible

**Interpretation:** a simulation can list bands that a real gel cannot resolve.
Small fragments may run off the gel or stain weakly.

**Actions:** redesign the amplicon, adapt agarose percentage and ladder, or
consider another method. Changing labels in a virtual gel cannot overcome a
physical resolution limit.

## Tm differs from another tool

**Possible causes:** different Na+, Mg2+, dNTP, or oligo concentrations,
another nearest-neighbour table, or another salt correction.

**Actions:** match all conditions before comparing. Use
`compare_tm_conditions()` to quantify the difference and validate with the real
PCR chemistry.

## ARMS-PCR has no valid sets or ambiguous bands

Inspect `design$diagnostics` and `design$exclusion_diagnostics`; confirm alleles
and orientation; and adjust windows or sizes only within justified limits. In
the laboratory, optimise temperature, primer ratio, and the deliberate
mismatch. See the [ARMS-PCR guide](ARMS_PCR_USER_GUIDE.md); a virtual gel does
not guarantee allelic discrimination.

## An experimental pattern disagrees with the simulation

Review controls, DNA identity and quality, digestion, contamination, agarose
settings, and lane annotations. Repeat according to a predefined criterion and
confirm relevant results independently. The simulation is a reference for
investigating a discrepancy, not a reason to override experimental controls.
