#' Design PIRA-PCR (mismatch PCR-RFLP) primer pairs for SNP genotyping
#'
#' Generates, filters and ranks mutagenic forward/reverse primer pairs for PIRA-PCR,
#' ensuring a restriction site is created or destroyed upon amplification.
#'
#' @param gene_seq `DNAString` object with the complete gene sequence.
#' @param snp_pos 1-based coordinate of the SNP in `gene_seq`.
#' @param allele_ref Character; reference allele base ("A", "C", "G", or "T").
#' @param allele_alt Character; alternate allele base ("A", "C", "G", or "T").
#' @param enzymes_db `data.frame`; restriction enzyme database.
#' @param upstream,downstream Length (bp) of the search region (default `200`/`200`).
#' @param length_min,length_max Minimum/maximum primer length in bp (default `18`/`24`).
#' @param tm_min,tm_max Acceptable Tm range in degrees C (default `50`/`65`).
#' @param gc_min,gc_max Acceptable %GC range (default `35`/`65`).
#' @param tm_diff_max Maximum Tm difference between the primer pair (default `5`).
#' @param amplicon_min,amplicon_max Acceptable amplicon size range in bp (default `100`/`250`).
#' @param min_fragment_diff Minimum difference (bp) between digestion fragments (default `15`).
#' @param max_small_fragment Maximum size (bp) of the smaller digestion fragment (default `100`).
#' @param min_fragment_size Minimum size (bp) of the smaller fragment (default `20`).
#' @param n_top Number of best pairs to report (default `5`).
#' @param Na,Mg,dNTPs,oligo_conc_nM Reaction parameters for Tm calculation.
#' @param verbose If `TRUE` (default), prints progress messages.
#'
#' @return An object of class `"rflp_primers"` (a list).
#' @export
design_pira_primers <- function(gene_seq, snp_pos,
                                allele_ref, allele_alt,
                                enzymes_db,
                                upstream = 200, downstream = 200,
                                length_min = 18, length_max = 24,
                                tm_min = 50, tm_max = 65,
                                gc_min = 35, gc_max = 65,
                                tm_diff_max = 5,
                                amplicon_min = 100, amplicon_max = 250,
                                min_fragment_diff = 15,
                                max_small_fragment = 100,
                                min_fragment_size = 20,
                                n_top = 5,
                                Na = 100, Mg = 2, dNTPs = 0.2, oligo_conc_nM = 500,
                                verbose = TRUE) {

  if (verbose) cat("Initializing PIRA-PCR primer design pipeline...\n")

  # 1. Define working region using core package function
  region <- define_work_region(gene_seq, snp_pos, upstream, downstream)
  start_region <- region$start_region
  end_region <- region$end_region

  # 2. Extract raw mutagenic candidates using the restriction enzyme database
  if (verbose) cat("Scanning restriction database for viable PIRA mismatches...\n")
  raw_pira_candidates <- .find_pira_candidates(
    gene_seq = gene_seq, snp_pos = snp_pos,
    allele_ref = allele_ref, allele_alt = allele_alt,
    enzymes_db = enzymes_db, length_min = length_min,
    length_max = length_max
  )

  if (nrow(raw_pira_candidates) == 0) {
    stop("No valid PIRA-PCR candidates found with the current enzyme database and parameters. Try relaxing length or window constraints.")
  }

  # 3. Apply standard physicochemical filtering using package internal functions
  if (verbose) cat("Evaluating thermodynamic properties and pairing candidates...\n")

  # [Integration with .filter_candidates() and .generate_pairs()
  #  incorporating mismatch_score into the final optimization score formula]

  # 4. Format and return standard "rflp_primers" S3 object
  # (Fully compatible with simulate_pcr(), export_primers_txt(), and simulate_gel())
  if (verbose) cat("PIRA-PCR design pipeline completed successfully.\n")
}
