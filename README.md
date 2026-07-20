# `rflpSNP` User Guide

### In silico design and simulation of PCR-RFLP assays for SNP genotyping

This guide is written for someone who **has never used R**. If you already
know R, feel free to jump straight to section 4 (workflow).

---

## 1. What does `rflpSNP` do?

A **PCR-RFLP** (*Restriction Fragment Length Polymorphism*) assay lets you
determine a person's genotype for a known SNP without sequencing: you PCR-
amplify the region of the gene containing the SNP, cut the amplicon with a
restriction enzyme whose recognition site coincides with the SNP position,
and observe the resulting band pattern on an agarose gel. If the SNP
creates or destroys the cut site, each allele produces a different
pattern:

| Genotype | Expected gel pattern |
|---|---|
| Homozygous without the cut site | 1 band (full, uncut amplicon) |
| Heterozygous | 3 bands (full amplicon + the 2 fragments) |
| Homozygous with the cut site | 2 bands (the 2 fragments) |

`rflpSNP` automates **the entire design of this assay**, from the
reference sequence to the expected gel image, without needing commercial
software (Primer3Plus, SnapGene, etc.). The package does not replace lab
work: it gives you a candidate design that you must then confirm
experimentally (and, for the primers, verify with a tool such as IDT
OligoAnalyzer before ordering synthesis).

---

## 2. Installation

### 2.1 Requirements

- **R** >= 4.1, and preferably **RStudio** (makes everything below
  easier).
- The packages `rflpSNP` depends on some other R packages that comes from **CRAN**
  (`TmCalculator`, `ggplot2`) and others from **Bioconductor**
  (`Biostrings`, `S4Vectors`) — Bioconductor is a repository separate from
  CRAN and needs its own installer.

### 2.2 Install the dependencies

Paste this into the R console **once**:

```r
# CRAN dependencies
install.packages(c("TmCalculator", "ggplot2", "devtools"))

# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Biostrings", "S4Vectors"))
```

### 2.3 Install `rflpSNP`

Install it with:

```r
devtools::install_github("Mike-EA/rflpSNP")
```

### 2.4 Load the package

At the start of every working session:

```r
library(rflpSNP)
```

### 2.5 Verify that everything works

`TmCalculator`'s interface has changed across versions, so `rflpSNP`
includes a diagnostic function. Run it once after installing:

```r
check_tm_backend()
```

If you see `[OK] calc_tm() is working correctly; test Tm = ...`, you're
all set. If you see a `[DIAGNOSTIC]` message, copy that output — it means
your installed version of `TmCalculator` renamed some argument and
`calc_tm()` (file `R/tm_utils.R`) needs to be adjusted accordingly.

---

## 3. Before writing code: prepare your 3 inputs

`rflpSNP` needs three things that **you** must obtain for your gene of
interest:

### 3.1 The reference sequence (FASTA file)

A `.fa`/`.fasta` text file with the gene (or genomic region) sequence
containing the SNP. It looks like this:

```
>NC_000001.11 MTHFR gene region
ATGGTGTCTGCGGGAGTCTGCAGTTCCCGGTGTAAAATCAGGGCAGTGACGGCAGCTGT...
```

You can download it, for example, from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore)
by searching for the gene and using "Send to > File > FASTA".

### 3.2 The dbSNP flanking sequence

This is the text immediately **preceding** the SNP on the reference
strand. It comes from the SNP's record page on [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
(NIH); search for the rs ID (e.g. `rs1801133`) and copy the 5' flanking
sequence. `rflpSNP` uses it as an "anchor" to find the exact SNP
coordinate within your FASTA, without you having to count bases by hand.

### 3.3 The restriction enzyme motif

The recognition site of the enzyme you plan to use, in IUPAC code (e.g.
`GANTC` for *HinfI*, where `N` = any base). You'll find it in the enzyme's
technical data sheet (NEB, Thermo Fisher) or on
[REBASE](http://rebase.neb.com/rebase/rebase.html). You also need to know
**exactly where it cuts** within the motif (for *HinfI*, it cuts between
the 1st and 2nd base: `G^ANTC`).

> **Useful IUPAC codes:** `N` = any base, `R` = A or G, `Y` = C or T,
> `W` = A or T, `S` = C or G. `matchPattern()` (used internally by
> `rflpSNP`) interprets these automatically.

---

## 4. Complete step-by-step workflow

We'll use the case that inspired this package as our example: the
**MTHFR** gene and the **C677T** SNP (rs1801133), cut with the **HinfI**
enzyme.

### Step 1 — Load the reference sequence

```r
gene_seq <- read_gene_fasta("MTHFR_completeseq.fa")
```

Expected output:
```
Sequence loaded. Length: XXXX bp
```

`gene_seq` is now a `DNAString` object (a special Bioconductor class for
DNA sequences, much more efficient than a plain `character` for this kind
of analysis).

### Step 2 — Locate the SNP with the flanking sequence

```r
snp <- locate_snp(gene_seq, flank_seq = "GAAAAGCTGCGTGATGATGAAATCG")
```

Expected output:
```
Nucleotide at position 9644 (C) [found on the forward strand]
```

`snp` is a list with three elements:

```r
snp$snp_pos    # SNP coordinate in gene_seq, e.g. 9644
snp$snp_base   # base found on the reference allele, e.g. "C"
snp$strand     # "forward" or "complementary"
```

Save `snp$snp_pos` — you'll use it in the next step.

> **What if I get "not found on either strand"?** Check that you copied
> the flanking sequence exactly as it appears in dbSNP (no spaces or line
> breaks) and that it corresponds to the same gene/assembly as your FASTA.

### Step 3 (optional, to understand the process) — Define the working region

`design_primers()` does this internally, but you can run it separately to
see what region will be explored:

```r
region <- define_work_region(gene_seq, snp_pos = snp$snp_pos,
                              upstream = 160, downstream = 200)
```

`upstream`/`downstream` are **intentionally different** (160 vs. 200):
this keeps the SNP from landing exactly in the middle of the working
region, which later helps the two digestion fragments end up at
distinguishable sizes on the gel.

### Step 4 — Design the primers

This is the package's main function:

```r
design <- design_primers(gene_seq, snp_pos = snp$snp_pos)
```

With the default values, `design_primers()` searches for 18-24 bp
primers, Tm 50-65 degrees C, %GC 35-65%, with a 150-300 bp amplicon, and
requires the two digestion fragments to differ by at least 25 bp and the
small fragment to be between 40 and 100 bp (so it separates well on the
gel). All of these values are adjustable — see the table in section 4.1.

Expected output (abridged):
```
Working region: nucleotide 9484 to 9844 (length: 361 bp; SNP at relative position 161)
Generating forward and reverse primer candidates...
Raw forward candidates: 2408 | Raw reverse candidates: 2408
Forward candidates passing physicochemical filtering: 34
Reverse candidates passing physicochemical filtering: 29

--- PAIR FILTERING SUMMARY ---
Failed due to amplicon size:        412
Failed due to Tm difference:        58
Failed due to (hetero)dimer:        11
Failed due to RFLP gel resolution:  25
Successful pairs:                   9

=== Best primer pair selected ===
Forward: TGGTCTCTTCATCCCTCGCCTTGAA (Tm=58.9C, GC=48.0%)
Reverse: GTCAGCCTCAAAGAAAAGCTGCGTG (Tm=59.2C, GC=48.0%)
Expected amplicon: 184 bp (positions 9542-9725) | Estimated RFLP fragments: 146 + 38 bp
```

`design` is an object that stores everything you need:

```r
design                 # prints a summary table of the 5 best pairs
design$top             # full data.frame with the 5 best pairs
design$best_pair       # the best pair (first row of design$top)
design$n_pairs_total   # how many pairs in total met every criterion
```

Most important columns of `design$top`:

| Column | Meaning |
|---|---|
| `forward_seq`, `reverse_seq` | Primer sequence, 5' -> 3' |
| `forward_tm`, `reverse_tm` | Melting temperature (degrees C) |
| `forward_gc`, `reverse_gc` | GC content (%) |
| `amplicon_bp` | Expected amplicon size (bp) |
| `tm_diff` | Tm difference between the two primers |
| `heterodimer_dg` | Heuristic heterodimer risk (closer to 0 = better) |
| `fragment_1`, `fragment_2` | Estimated fragment sizes after the cut |
| `score` | Internal ranking index (lower = better) |
| `recommended` | `"YES"` only on the first row |

> **What do I do if `design_primers()` errors with "No valid pairs were
> found"?** Relax the most restrictive criteria first: `amplicon_max`,
> `tm_diff_max`, `min_fragment_diff` or `max_small_fragment`. See the
> parameter table below.

#### 4.1 `design_primers()` parameters you'll probably want to adjust

| Parameter | Default | What it controls |
|---|---|---|
| `upstream`, `downstream` | `160`, `200` | Size of the search region on each side of the SNP |
| `length_min`, `length_max` | `18`, `24` | Primer length (bp) |
| `tm_min`, `tm_max` | `50`, `65` | Acceptable Tm range (degrees C) |
| `gc_min`, `gc_max` | `35`, `65` | Acceptable %GC range |
| `tm_diff_max` | `5` | Maximum Tm difference between forward and reverse |
| `amplicon_min`, `amplicon_max` | `150`, `300` | Acceptable amplicon size (bp) |
| `min_fragment_diff` | `25` | Minimum difference (bp) between the two cut fragments |
| `max_small_fragment` | `100` | Maximum size the small fragment may have |
| `min_fragment_size` | `40` | Minimum size of the small fragment (so it's visible on a gel) |
| `n_top` | `5` | How many of the best pairs to report |

Example relaxing criteria because the default design found no pairs:

```r
design <- design_primers(gene_seq, snp_pos = snp$snp_pos,
                          amplicon_max = 400,
                          min_fragment_diff = 15)
```

#### 4.2 If the Tm looks different from SnapGene / IDT OligoAnalyzer

This is expected — see the "Details" section of `?calc_tm` for the full
explanation (mainly: whether Mg2+ is included in the salt conditions).
`compare_tm_conditions()` lets you see, in one call, how much the Tm of
the same primer shifts across a few typical condition profiles:

```r
compare_tm_conditions(design$best_pair$forward_seq)
```

### Step 5 — Export the report to text

```r
export_primers_txt(design, output_file = "primers_MTHFR_C677T.txt")
```

This creates a `.txt` file with the best pair in detail and the summary
table of the 5 best candidates — handy for pasting into your lab notebook
or sharing with students.

### Step 6 — Simulate the PCR in silico

Using the best primer pair, we generate the real amplicon (cut out of the
reference sequence):

```r
bp <- design$best_pair

pcr <- simulate_pcr(gene_seq,
                     fwd_primer = bp$forward_seq,
                     rev_primer = bp$reverse_seq)
```

Expected output:
```
Amplicon: 184 bp (positions 9542-9725)
Binding sites (<=2 mismatches) - Fwd: 1 | Rev: 1
```

If you see `Fwd: 2` or more, there is more than one possible binding site
in your reference sequence — check that primer's specificity (e.g. with a
BLAST against the full genome) before continuing.

`pcr` stores:
```r
pcr$amplicon   # the amplicon as a DNAString
pcr$start      # start coordinate in gene_seq
pcr$end        # end coordinate in gene_seq
pcr$size       # size in bp
```

### Step 7 — Locate the restriction enzyme site in the amplicon

```r
site <- find_restriction_site(pcr, enzyme_motif = "GANTC", cut_offset = 1)
```

Expected output:
```
Found 1 site(s) for 'GANTC' in the amplicon.
```

`cut_offset = 1` means *HinfI* cuts right after the first base of the
motif (`G^ANTC`). If you use a different enzyme, adjust this value
according to where that specific enzyme cuts.

> If `find_restriction_site()` reports **0 sites**, that's not necessarily
> an error: it can mean that, on the allele you loaded in your FASTA (the
> reference allele), the cut site simply **doesn't exist** — i.e. the SNP
> is exactly what *creates* the site on the alternate allele (see Step 9).

### Step 8 — Plot the amplicon and the sequence

Schematic map (primers + cut site):

```r
plot_amplicon_map(pcr, bp$forward_seq, bp$reverse_seq, site,
                   enzyme_name = "HinfI")
```

Letter-by-letter view of the full amplicon, with the site highlighted:

```r
plot_sequence_map(pcr, site, enzyme_name = "HinfI")
```

Both functions return a `ggplot` object, so you can save them as images:

```r
map <- plot_amplicon_map(pcr, bp$forward_seq, bp$reverse_seq, site, enzyme_name = "HinfI")
ggplot2::ggsave("amplicon_map.png", plot = map, width = 8, height = 4, dpi = 300)
```

### Step 9 — Get the fragments of the **other** allele (for the full gel)

A PCR-RFLP assay compares **two alleles**: the reference one (the one in
your FASTA) and the alternate one (the one that creates or removes the cut
site). `rflpSNP` always works on whatever sequence you give it, so to
simulate the other allele you substitute the SNP base and repeat the
simulation:

```r
# Copy the sequence and change the SNP base to the alternate allele (e.g. "T")
gene_seq_alt <- gene_seq
Biostrings::subseq(gene_seq_alt,
                    start = snp$snp_pos, end = snp$snp_pos) <- Biostrings::DNAString("T")

pcr_alt <- simulate_pcr(gene_seq_alt, bp$forward_seq, bp$reverse_seq)
site_alt <- find_restriction_site(pcr_alt, enzyme_motif = "GANTC", cut_offset = 1)
```

One of the two alleles (reference or alternate) will have the cut site and
the other won't — compare `site$n_sites` against `site_alt$n_sites` to
confirm which is which, and note down the fragment sizes (`site$sites`) of
the allele that does get cut.

### Step 10 — Simulate the agarose gel

Using the amplicon size and the two fragments from the allele that does
get cut:

```r
simulate_gel(amplicon_size = pcr$size,
             fragment_sizes = c(146, 38),
             genotype_labels = c("Homozygous C/C", "Heterozygous C/T", "Homozygous T/T"))
```

This produces the virtual gel image with 4 lanes: molecular weight marker,
uncut homozygote, heterozygote (3 bands) and cut homozygote (2 bands).
Save it just like the other plots with `ggplot2::ggsave()`.

### Step 11 (shortcut) — The whole workflow in one call

For a quick first look (e.g. so a student can see the full workflow right
away), `run_pcr_rflp_pipeline()` chains steps 1(b) through 7
automatically:

```r
result <- run_pcr_rflp_pipeline(
  gene_seq, flank_seq = "GAAAAGCTGCGTGATGATGAAATCG",
  enzyme_motif = "GANTC", cut_offset = 1
)

result$snp
result$design
result$pcr_result
result$restriction_result
```

For classroom use, **running each step separately the first time**
(previous sections) is recommended, and using this shortcut only once the
workflow is understood.

---

## 5. Quick glossary for bioinformatics beginners

- **SNP** (*Single Nucleotide Polymorphism*): a genome position where
  different people can have different bases (e.g. C or T).
- **Primer**: short oligonucleotide (typically 18-24 bp) that serves as the
  starting point for the polymerase to copy DNA during PCR.
- **Tm** (melting temperature): the temperature at which half of a
  primer's molecules are bound to their complementary sequence and the
  other half are separated. Two primers of the same pair should have
  similar Tm so both bind well at the same PCR temperature.
- **%GC**: percentage of G and C bases in the primer; affects binding
  stability.
- **Amplicon**: the DNA fragment resulting from PCR (the region between
  the forward and reverse primer, both included).
- **Dimer**: unwanted pairing between two copies of the same primer
  (self-dimer) or between the forward and reverse primer (heterodimer),
  which consumes reagents and reduces PCR yield.
- **Hairpin**: a primer that folds back on itself because part of its
  sequence is complementary to another part of the same sequence.
- **Restriction enzyme**: a protein that cuts DNA at a specific sequence
  site (e.g. *HinfI* cuts at `GANTC`).
- **RFLP**: difference in restriction fragment size between individuals,
  caused by a SNP that creates or destroys a cut site.

---

## 6. Best practices before going to the lab

1. **Verify the primers with a specialized tool** (e.g. IDT OligoAnalyzer)
   before ordering synthesis: `rflpSNP`'s Tm, dimer and hairpin estimates
   are fast heuristics for *comparing candidates against each other*, not
   definitive thermodynamic calculations.
2. **Confirm primer specificity** (e.g. with a BLAST against the full
   genome/transcriptome), especially if `simulate_pcr()` reported more
   than one potential binding site.
3. **Verify the enzyme's cut site against the supplier's technical
   sheet**; `cut_offset` depends on each enzyme.

---

## 7. Where to find each function (package map)

| File in `R/` | Functions | What it does |
|---|---|---|
| `read_gene_fasta.R` | `read_gene_fasta()` | Loads the reference sequence |
| `locate_snp.R` | `locate_snp()` | Locates the SNP using the flanking sequence |
| `define_work_region.R` | `define_work_region()` | Defines the primer search region |
| `tm_utils.R` | `calc_tm()`, `check_tm_backend()`, `compare_tm_conditions()` | Tm calculation, `TmCalculator` diagnostics, and condition comparison |
| `dimer_hairpin.R` | `evaluate_dimer()`, `evaluate_hairpin()` | Dimer and hairpin risk |
| `primer_candidates.R` | *(internal)* | Candidate/pair generation and filtering |
| `design_primers.R` | `design_primers()` | Complete primer design |
| `export_primers.R` | `export_primers_txt()` | Exports the report to `.txt` |
| `pcr_insilico.R` | `simulate_pcr()` | Simulates the PCR and returns the amplicon |
| `restriction_site.R` | `find_restriction_site()` | Locates the enzyme site in the amplicon |
| `plot_amplicon.R` | `plot_amplicon_map()`, `plot_sequence_map()` | Amplicon plots |
| `gel_simulation.R` | `simulate_gel()` | Virtual agarose gel |
| `pipeline.R` | `run_pcr_rflp_pipeline()` | Shortcut for the whole workflow |

---

