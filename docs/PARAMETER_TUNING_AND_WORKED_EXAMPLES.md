# Parameter tuning and worked examples

Use parameters to explore a design space, not to convert a weak candidate into
an experimentally valid assay. Record every non-default setting with the FASTA,
assembly, alleles, package version, and exported report.

## Parameters by workflow stage

| Functions | Key parameters | Effect |
|---|---|---|
| `read_gene_fasta()`, `locate_snp()` | `seq_index`, `flank_seq`, `snp_offset` | Select a correct reference and locate the SNP. Use a longer unique flank before changing `snp_offset`. |
| `calc_tm()`, `compare_tm_conditions()` | `Na`, `Mg`, `dNTPs`, `oligo_conc_nM`, `nn_table`, `salt_corr_method`, `profiles` | Match Tm calculations to the intended PCR chemistry. |
| `define_work_region()` | `upstream`, `downstream`, `max_region_size` | Defines the search geometry. A wider or asymmetric window can create options but does not repair poor specificity. |
| `design_primers()` | length, Tm, GC, and `tm_diff_max` limits | Filter individual PCR-RFLP primers. Relax one criterion at a time and document why. |
| `design_primers()` | amplicon, fragment, `max_candidates_per_strand`, `n_top` limits | Controls gel-resolvable products, runtime, and alternatives reported. A small cap is exploratory, not exhaustive. |
| `simulate_pcr()` | `max_mismatch` | Counts approximate sites in the supplied FASTA; it is not a genome-wide specificity search. |
| `find_restriction_site()` | `enzyme_motif`, `cut_offset` | Defines the recognition motif and cut coordinate; confirm both from enzyme documentation. |
| `simulate_gel()` | sizes, `ladder_sizes`, `genotype_labels` | Changes an idealised visual model only. |
| `design_arms_primers()` | `outer_flank`, outer/inner lengths, `mismatch_positions` | Controls four-primer geometry and intentional 3' discrimination. Mismatch choices require laboratory optimisation. |
| `design_arms_primers()` | Tm/GC chemistry, dimer/hairpin thresholds, product ranges, `min_band_diff` | Balances primer compatibility and gel readability. Relaxed structure thresholds require independent review. |
| `design_arms_primers()` | candidate-pool caps and `n_top` | Bounds exploration time; repeat with broader pools before a final recommendation. |
| `plot_*()`, `export_*()` | labels, site index, row width, output path | Improve interpretation and preserve traceability; they do not change the assay. |
| `run_*_pipeline()` | `export_txt`, `output_file`, `...` | Runs a workflow and forwards `...` to the design function. Use a stepwise workflow first when tuning. |

`evaluate_dimer()` and `evaluate_hairpin()` are inspection tools. Their
heuristic results prioritise candidates and do not replace specialised
thermodynamic software.

## MTHFR C677T PCR-RFLP example

The documented MTHFR assay has a 323 bp amplicon and HinfI fragments of 225 and
98 bp. These values favour bands that can be resolved with an appropriate gel.

The following reproducible workflow uses the bundled GRCh38 record, selects a
323 bp amplicon, and verifies the actual HinfI cut after modelling the C→T
allelic change:

```r
library(rflpSNP)

# The file contains two assemblies. Select the documented GRCh38 record.
reference <- read_gene_fasta(
  system.file("extdata", "MTHFR_complete.fna", package = "rflpSNP"),
  seq_index = 1
)
snp <- locate_snp(reference, "GAAAAGCTGCGTGATGATGAAATCG")
stopifnot(snp$snp_pos == 9644L, snp$snp_base == "C")

# A 300--400 bp product and strict fragment constraints favour gel-resolvable
# products. A cap of 100 candidates per strand keeps the final search broad.
design <- design_primers(
  reference, snp_pos = snp$snp_pos,
  upstream = 140L, downstream = 300L, min_distance_to_snp = 20L,
  length_min = 18L, length_max = 20L,
  tm_min = 50, tm_max = 65, gc_min = 35, gc_max = 65, tm_diff_max = 4,
  amplicon_min = 300L, amplicon_max = 400L,
  min_fragment_size = 80L, max_small_fragment = 140L,
  min_fragment_diff = 100L, max_candidates_per_strand = 100L, n_top = 10L,
  Na = 100, Mg = 2, dNTPs = 0.2, oligo_conc_nM = 500, verbose = TRUE
)

best_pair <- design$best_pair
pcr <- simulate_pcr(reference, best_pair$forward_seq, best_pair$reverse_seq)
stopifnot(pcr$size == best_pair$amplicon_bp)

# Construct the alternate T allele only at the SNP and find the real cut.
snp_in_amplicon <- snp$snp_pos - pcr$start + 1L
alternate_amplicon <- pcr$amplicon
Biostrings::subseq(alternate_amplicon, snp_in_amplicon, snp_in_amplicon) <-
  Biostrings::DNAString("T")

reference_site <- find_restriction_site(pcr, "GANTC", cut_offset = 1)
alternate_site <- find_restriction_site(alternate_amplicon, "GANTC", 1)
stopifnot(reference_site$n_sites == 0L, alternate_site$n_sites == 1L)

cut_position <- alternate_site$sites$cut_pos[1]
fragments <- sort(c(cut_position, pcr$size - cut_position), decreasing = TRUE)
stopifnot(identical(fragments, c(225, 98)))

gel <- simulate_gel(pcr$size, fragments, genotype_labels = c("C/C", "C/T", "T/T"))
gel
```

![MTHFR C677T virtual gel](images/mthfr-c677t-virtual-gel.png)

The recommended pair is `TGACTGTCATCCCTATTG` / `GAAGAACTCAGCGAACTC`. The
heterozygote has 323, 225, and 98 bp bands. This is an output from
`simulate_gel()`, not an experimental gel. Confirm enzyme conditions and
primer specificity before synthesis.

## FTO rs8050136 tetra-primer ARMS-PCR example

The bundled FTO regression case uses a 171 bp reference window, SNP coordinate
86, C reference allele, A alternate allele, 20 bp outer primers, 18--20 bp
inner primers, Tm 55--72 °C, and a minimum 20 bp band difference. Its expected
products are 152 bp control, 105 bp reference, and 84 bp alternate.

The following complete code generated the figure. It uses the fixed reference
window and all exploratory parameters from the regression case:

```r
library(rflpSNP)

# Load the bundled full FTO reference and locate rs8050136 reproducibly.
reference <- read_gene_fasta(
  system.file("extdata", "FTO_complete.fa", package = "rflpSNP")
)
snp <- locate_snp(reference, "CATGCCAGTTGCCCACTGTGGCAAT")
stopifnot(snp$snp_pos == 78401L, snp$snp_base == "C")

# Tight product ranges and a 20 bp minimum difference make the three bands
# readable; capped pools keep this example quick to reproduce.
design <- design_arms_primers(
  reference, snp_pos = snp$snp_pos, ref_allele = snp$snp_base, alt_allele = "A",
  outer_flank = 85L, outer_length_min = 20L, outer_length_max = 20L,
  inner_length_min = 18L, inner_length_max = 20L,
  tm_min = 55, tm_max = 72, gc_min = 35, gc_max = 65,
  dimer_dg_min = -8, hairpin_dg_min = -6, heterodimer_dg_min = -9,
  control_amplicon_min = 150L, control_amplicon_max = 200L,
  allele_amplicon_min = 80L, allele_amplicon_max = 150L,
  min_band_diff = 20L, max_candidates_per_pool = 8L,
  max_raw_candidates_per_pool = 8L
)

# Confirm the three predicted genotypes and render the virtual gel.
pcr_result <- simulate_arms_pcr(
  reference, design$best_set, snp_pos = snp$snp_pos,
  ref_allele = snp$snp_base, alt_allele = "A"
)

pcr_result$products
simulate_arms_gel(pcr_result)
```

![FTO rs8050136 virtual gel](images/fto-rs8050136-arms-virtual-gel.png)

The heterozygote shows all three bands. The regression thresholds are an
exploratory example and require independent review before experimental use.
