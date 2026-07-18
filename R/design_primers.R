#' Design primer pairs for PCR-RFLP
#'
#' Generates, filters and ranks forward/reverse primer pairs around a SNP,
#' taking into account standard physicochemical properties (Tm, %GC), risk
#' of self-dimer, heterodimer and hairpin formation, amplicon size, and the
#' expected resolution of the RFLP fragments on an agarose gel (so that the
#' SNP does not end up centered in the amplicon and the two cut fragments
#' are distinguishable).
#'
#' @param gene_seq `DNAString` object with the complete gene sequence.
#' @param snp_pos 1-based coordinate of the SNP in `gene_seq` (see
#'   [locate_snp()]).
#' @param upstream,downstream Length (bp) of the primer search region
#'   upstream/downstream of the SNP (default `160`/`200`; see
#'   [define_work_region()]).
#' @param min_distance_to_snp Minimum distance (bp) between the SNP and the
#'   nearest edge of each primer, to leave room for the restriction
#'   enzyme's recognition site (default `20`).
#' @param length_min,length_max Minimum/maximum primer length in bp (default
#'   `18`/`24`).
#' @param tm_min,tm_max Acceptable Tm range in degrees C (default
#'   `50`/`65`).
#' @param gc_min,gc_max Acceptable %GC range (default `35`/`65`).
#' @param tm_diff_max Maximum Tm difference between the forward and reverse
#'   primer of the same pair (default `5`).
#' @param amplicon_min,amplicon_max Acceptable amplicon size range in bp
#'   (default `150`/`300`, typical for PCR-RFLP).
#' @param min_fragment_diff Minimum difference (bp) between the two
#'   restriction fragments for them to be considered distinguishable on a
#'   gel (default `25`).
#' @param max_small_fragment Maximum size (bp) the smaller digestion
#'   fragment may have for the pair to be considered suitable (default
#'   `100`; small fragments separate better visually on the logarithmic
#'   migration scale of an agarose gel).
#' @param min_fragment_size Minimum size (bp) of the smaller fragment, to
#'   avoid fragments too small to visualize on a standard gel (default
#'   `40`).
#' @param max_candidates_per_strand Cap on preselected candidates per strand
#'   before generating all fwd x rev combinations, to bound execution time
#'   (default `40`).
#' @param n_top Number of best pairs to report (default `5`).
#' @param Na,Mg,dNTPs,oligo_conc_nM PCR reaction parameters used for the Tm
#'   calculation (see [calc_tm()]).
#' @param verbose If `TRUE` (default), prints progress messages and
#'   filtering summaries at each stage.
#'
#' @return An object of class `"rflp_primers"` (a list) with:
#'   \item{top}{`data.frame` with the `n_top` best pairs, including a
#'     `recommended` column.}
#'   \item{best_pair}{The best pair (first row of `top`).}
#'   \item{n_candidates_forward, n_candidates_reverse}{Number of individual
#'     candidates that passed the physicochemical filtering.}
#'   \item{n_pairs_total}{Total number of pairs that met all criteria.}
#'   \item{parameters}{List with the parameters used, for traceability.}
#'
#' @seealso [export_primers_txt()] to export the results to a text file,
#'   [define_work_region()] for the earlier step of delimiting the search
#'   region.
#'
#' @export
design_primers <- function(gene_seq, snp_pos,
                            upstream = 160, downstream = 200,
                            min_distance_to_snp = 20,
                            length_min = 18, length_max = 24,
                            tm_min = 50, tm_max = 65,
                            gc_min = 35, gc_max = 65,
                            tm_diff_max = 5,
                            amplicon_min = 150, amplicon_max = 300,
                            min_fragment_diff = 25,
                            max_small_fragment = 100,
                            min_fragment_size = 40,
                            max_candidates_per_strand = 40,
                            n_top = 5,
                            Na = 100, Mg = 2, dNTPs = 0.2, oligo_conc_nM = 500,
                            verbose = TRUE) {

  region <- define_work_region(gene_seq, snp_pos, upstream, downstream)
  start_region <- region$start_region
  end_region <- region$end_region

  limit_forward <- snp_pos - min_distance_to_snp
  if (limit_forward < start_region) {
    stop("The region upstream of the SNP is too short to design a forward primer. Increase 'upstream' or reduce 'min_distance_to_snp'.")
  }
  fwd_search_region <- Biostrings::subseq(gene_seq, start_region, limit_forward)

  start_reverse <- snp_pos + min_distance_to_snp
  if (start_reverse > end_region) {
    stop("The region downstream of the SNP is too short to design a reverse primer. Increase 'downstream' or reduce 'min_distance_to_snp'.")
  }
  rev_search_region <- Biostrings::subseq(gene_seq, start_reverse, end_region)

  if (verbose) cat("Generating forward and reverse primer candidates...\n")
  candidates_forward_raw <- .extract_candidates(fwd_search_region, length_min, length_max)
  candidates_reverse_raw <- .extract_candidates(rev_search_region, length_min, length_max)

  if (verbose) {
    cat(sprintf(
      "Raw forward candidates: %d | Raw reverse candidates: %d\n",
      length(candidates_forward_raw), length(candidates_reverse_raw)
    ))
  }

  candidates_f <- .filter_candidates(
    candidates_forward_raw, region_offset = start_region, is_reverse = FALSE,
    tm_min = tm_min, tm_max = tm_max, gc_min = gc_min, gc_max = gc_max,
    Na = Na, Mg = Mg, dNTPs = dNTPs, oligo_conc_nM = oligo_conc_nM
  )
  candidates_r <- .filter_candidates(
    candidates_reverse_raw, region_offset = start_reverse, is_reverse = TRUE,
    tm_min = tm_min, tm_max = tm_max, gc_min = gc_min, gc_max = gc_max,
    Na = Na, Mg = Mg, dNTPs = dNTPs, oligo_conc_nM = oligo_conc_nM
  )

  if (verbose) {
    cat(sprintf("Forward candidates passing physicochemical filtering: %d\n", nrow(candidates_f)))
    cat(sprintf("Reverse candidates passing physicochemical filtering: %d\n", nrow(candidates_r)))
  }

  if (nrow(candidates_f) == 0 || nrow(candidates_r) == 0) {
    stop("No valid candidates were found for one of the two strands. Relax the criteria (Tm, GC, dimers, hairpins).")
  }

  candidates_f <- .preselect_candidates(candidates_f, max_candidates_per_strand, tm_min, tm_max, gc_min, gc_max)
  candidates_r <- .preselect_candidates(candidates_r, max_candidates_per_strand, tm_min, tm_max, gc_min, gc_max)

  candidate_pairs <- .generate_pairs(
    candidates_f, candidates_r, snp_pos = snp_pos,
    amplicon_min = amplicon_min, amplicon_max = amplicon_max,
    tm_diff_max = tm_diff_max,
    min_fragment_diff = min_fragment_diff,
    max_small_fragment = max_small_fragment,
    min_fragment_size = min_fragment_size,
    verbose = verbose
  )

  if (nrow(candidate_pairs) == 0) {
    stop("No valid pairs were found. Consider relaxing amplicon_min/max, tm_diff_max, min_fragment_diff or max_small_fragment.")
  }

  tm_target <- mean(c(tm_min, tm_max))
  gc_target <- mean(c(gc_min, gc_max))

  candidate_pairs$score <- with(candidate_pairs,
    abs(forward_tm - tm_target) + abs(reverse_tm - tm_target) +
      abs(forward_gc - gc_target) * 0.2 + abs(reverse_gc - gc_target) * 0.2 +
      tm_diff * 1.5 + pmax(0, -heterodimer_dg - 4) * 1.0
  )

  sorted_pairs <- candidate_pairs[order(candidate_pairs$score), ]
  top <- utils::head(sorted_pairs, n_top)
  top$recommended <- c("YES", rep("NO", nrow(top) - 1))
  rownames(top) <- NULL

  result <- list(
    top = top,
    best_pair = top[1, ],
    n_candidates_forward = nrow(candidates_f),
    n_candidates_reverse = nrow(candidates_r),
    n_pairs_total = nrow(candidate_pairs),
    parameters = list(
      snp_pos = snp_pos, upstream = upstream, downstream = downstream,
      min_distance_to_snp = min_distance_to_snp,
      length_min = length_min, length_max = length_max,
      tm_min = tm_min, tm_max = tm_max, gc_min = gc_min, gc_max = gc_max,
      tm_diff_max = tm_diff_max, amplicon_min = amplicon_min, amplicon_max = amplicon_max,
      min_fragment_diff = min_fragment_diff,
      max_small_fragment = max_small_fragment,
      min_fragment_size = min_fragment_size
    )
  )
  class(result) <- "rflp_primers"

  if (verbose) {
    cat("\n=== Best primer pair selected ===\n")
    bp <- result$best_pair
    cat(sprintf("Forward: %s (Tm=%.1f\u00B0C, GC=%.1f%%)\n", bp$forward_seq, bp$forward_tm, bp$forward_gc))
    cat(sprintf("Reverse: %s (Tm=%.1f\u00B0C, GC=%.1f%%)\n", bp$reverse_seq, bp$reverse_tm, bp$reverse_gc))
    cat(sprintf(
      "Expected amplicon: %d bp (positions %d-%d) | Estimated RFLP fragments: %d + %d bp\n",
      bp$amplicon_bp, bp$forward_start, bp$reverse_end, bp$fragment_1, bp$fragment_2
    ))
  }

  result
}

#' @export
print.rflp_primers <- function(x, ...) {
  cat("<rflp_primers>\n")
  cat(sprintf(
    "%d total candidate pair(s) evaluated; showing the best %d.\n",
    x$n_pairs_total, nrow(x$top)
  ))
  print(x$top[, c(
    "forward_seq", "forward_tm", "reverse_seq", "reverse_tm",
    "amplicon_bp", "fragment_1", "fragment_2", "score", "recommended"
  )], row.names = FALSE)
  invisible(x)
}
