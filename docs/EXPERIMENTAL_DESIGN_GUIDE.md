# Experimental Design Guide for SNP Genotyping

This guide places `rflpSNP` computational analysis within experimental assay
design. The package generates reproducible candidates; biological judgement,
controls, and laboratory validation remain necessary for the final decision.

## 1. Define the question and reference

Record the SNP identifier, alleles to distinguish, organism, assembly, and
reference version before designing primers. Also record the FASTA identifier
and source. Coordinates can differ among assemblies, transcripts, and database
versions.

The FASTA should contain the relevant genomic region in one orientation. When
it contains several sequences, `read_gene_fasta()` uses the record selected by
`seq_index`; document that choice. The base returned by `locate_snp()` must
match the expected reference allele before continuing.

## 2. Choose a strategy

PCR-RFLP is appropriate when an allelic change creates or removes a restriction
site that yields distinguishable fragments. The target is amplified and
digested; in an ideal model, the uncut allele retains the full amplicon, the
cut allele yields two fragments, and a heterozygote shows both patterns.

Tetra-primer ARMS-PCR is an alternative when there is no useful restriction-site
change. It uses two outer primers for a control band and two inner primers whose
3' ends discriminate alleles. See the [ARMS-PCR user guide](ARMS_PCR_USER_GUIDE.md)
for its specific assumptions.

The absence of a valid design does not prove that an SNP cannot be genotyped.
It may mean that the criteria, available window, or selected strategy are not
appropriate. Record that decision rather than forcing a poor candidate.

## 3. Establish whether PCR-RFLP is viable

For every candidate enzyme, confirm with a primary source or supplier data
sheet:

- the recognition motif and cut positions on both strands;
- whether the SNP alters that motif for each allele;
- that the amplicon contains no additional confusing sites;
- methylation sensitivity, buffer requirements, temperature, and digestion time;
- that expected fragments can be resolved with the intended gel system.

`find_restriction_site()` searches an IUPAC motif and calculates a cut from the
supplied `cut_offset`. It does not identify an enzyme from its motif or verify
biochemical conditions.

## 4. Treat primer design as a working hypothesis

`design_primers()` searches both sides of the SNP and ranks candidates by
length, Tm, %GC, heuristic complementarity, amplicon size, and estimated RFLP
fragment sizes. Its parameters are exploration filters, not universal rules.
Change a constraint only for a documented reason:

- expand the window or amplicon size when the region offers few options;
- adapt the Tm range to the PCR chemistry and protocol being used;
- retain fragment sizes that can be visualised with the selected agarose and
  ladder;
- compare the leading candidates instead of relying only on the first row of
  `design$top`.

Tm depends on ion, dNTP, and oligo concentrations. `calc_tm()` uses a
nearest-neighbour model and the PCR conditions specified in the function. Use
`compare_tm_conditions()` to compare settings; other software may report a
different value when conditions are not identical.

## 5. Review specificity and products

`simulate_pcr()` builds an amplicon from the first exact matches of both primers
in the supplied sequence. It also counts matches with up to `max_mismatch`
differences as a local alert. This checks the input FASTA, not a genome,
transcriptome, pseudogene, or population-variant catalogue.

Before synthesis, evaluate each primer against the relevant genome or
transcriptome with an appropriate specificity tool. Review dimers, hairpins,
and 3' ends with specialised thermodynamic software. Package functions are
fast screens for comparing candidates.

## 6. Plan controls and interpretation

Use suitable controls: a no-template control, known genotypes when available,
a size ladder, and—for RFLP—a digestion control. Define beforehand which
patterns will be accepted, when an ambiguous band will be repeated, and how a
discrepancy will be confirmed.

Band intensity, partial digestion, degraded DNA, inhibitors, concentration
variation, and gel resolution can alter a real sample. The [limitations guide](LIMITATIONS_AND_VALIDATION.md)
explains these differences between simulation and laboratory work.

## 7. Preserve traceability

Store the FASTA, assembly version, flank, alleles, motif and cut position,
arguments, `rflpSNP` and `TmCalculator` versions, date, and exported reports.
The pipeline and export functions preserve useful results but do not replace an
experimental notebook.
