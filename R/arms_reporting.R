#' Export tetra-primer ARMS-PCR design results to text
#'
#' Writes the recommended four-primer set, the ranked candidates, predicted
#' genotype bands, parameters and experimental-use warnings.
#'
#' @param design A `rflp_arms_primers` object from [design_arms_primers()].
#' @param output_file Path to the report to create.
#'
#' @return Invisibly, the path of the written report.
#' @export
export_arms_primers_txt <- function(design, output_file = "arms_pcr_primers_results.txt") {
  if (!inherits(design, "rflp_arms_primers") || is.null(design$n_valid_sets)) {
    stop("'design' must be an object returned by design_arms_primers().")
  }
  if (!design$n_valid_sets) {
    stop("Cannot export an ARMS-PCR report because the design contains no valid set.")
  }
  best <- design$best_set
  report <- c(
    "===========================================================",
    "             TETRA-PRIMER ARMS-PCR DESIGN RESULTS",
    "===========================================================", "",
    "--- RECOMMENDED SET ---",
    sprintf("Layout: %s | score: %.3f", best$layout, best$score),
    sprintf("External control product: %d bp", best$control_amplicon_bp),
    sprintf("Reference-specific product: %d bp", best$ref_amplicon_bp),
    sprintf("Alternative-specific product: %d bp", best$alt_amplicon_bp), "",
    "Outer forward (5' -> 3'):",
    sprintf("  %s | FASTA %d-%d | Tm %.1f C | GC %.1f%%", best$outer_forward_seq,
            best$outer_forward_start, best$outer_forward_end, best$outer_forward_tm, best$outer_forward_gc),
    "Outer reverse (5' -> 3'):",
    sprintf("  %s | FASTA %d-%d | Tm %.1f C | GC %.1f%%", best$outer_reverse_seq,
            best$outer_reverse_start, best$outer_reverse_end, best$outer_reverse_tm, best$outer_reverse_gc),
    "Reference inner (5' -> 3'):",
    sprintf("  %s | %s | FASTA %d-%d | deliberate mismatch %s", best$ref_inner_seq,
            best$ref_inner_strand, best$ref_inner_start, best$ref_inner_end, best$ref_mismatch_3prime),
    "Alternative inner (5' -> 3'):",
    sprintf("  %s | %s | FASTA %d-%d | deliberate mismatch %s", best$alt_inner_seq,
            best$alt_inner_strand, best$alt_inner_start, best$alt_inner_end, best$alt_mismatch_3prime), "",
    "--- PREDICTED BANDS ---",
    sprintf("ref/ref: %d + %d bp (control + reference)", best$control_amplicon_bp, best$ref_amplicon_bp),
    sprintf("ref/alt: %d + %d + %d bp", best$control_amplicon_bp, best$ref_amplicon_bp, best$alt_amplicon_bp),
    sprintf("alt/alt: %d + %d bp (control + alternative)", best$control_amplicon_bp, best$alt_amplicon_bp), "",
    "--- TOP CANDIDATES ---"
  )
  summary <- utils::capture.output(print(design$top[, c(
    "layout", "control_amplicon_bp", "ref_amplicon_bp", "alt_amplicon_bp",
    "worst_cross_dimer_dg", "score", "recommended"
  )], row.names = FALSE, right = FALSE))
  parameters <- utils::capture.output(utils::str(design$parameters, give.attr = FALSE))
  report <- c(report, summary, "", "--- PARAMETERS USED ---", parameters, "",
    "--- IMPORTANT LIMITATIONS ---",
    "* This is an in-silico candidate design, not experimental validation.",
    "* Confirm genomic specificity and secondary structures before synthesis.",
    "* ARMS-PCR needs a polymerase without 3'->5' proofreading activity.",
    "* Optimize the deliberate mismatch and annealing conditions experimentally.")
  writeLines(report, con = output_file)
  message(sprintf("ARMS-PCR results exported to: %s", output_file))
  invisible(output_file)
}
