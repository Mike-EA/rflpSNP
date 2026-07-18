#' Export primer design results to a text file
#'
#' Writes a human-readable plain-text report with the recommended best
#' primer pair and a summary table of the `n_top` best candidates.
#'
#' @param design `"rflp_primers"` object returned by [design_primers()].
#' @param output_file Path to the output file (default
#'   `"primers_PCR_RFLP_results.txt"`).
#'
#' @return Invisibly, the path of the written file.
#' @export
export_primers_txt <- function(design, output_file = "primers_PCR_RFLP_results.txt") {
  if (!inherits(design, "rflp_primers")) {
    stop("'design' must be an object returned by design_primers().")
  }

  bp <- design$best_pair

  report <- c(
    "===========================================================",
    "            PCR-RFLP PRIMER DESIGN RESULTS                 ",
    "===========================================================\n",
    "--- BEST PAIR SELECTED (RECOMMENDED) ---",
    sprintf("Expected amplicon  : %d bp", bp$amplicon_bp),
    sprintf("Positions (FASTA)  : %d to %d", bp$forward_start, bp$reverse_end),
    sprintf("Estimated RFLP fragments : %d bp + %d bp", bp$fragment_1, bp$fragment_2),
    sprintf("Tm difference      : %.1f \u00B0C", bp$tm_diff),
    sprintf("Heterodimer risk (approx. Delta G): %.1f kcal/mol\n", bp$heterodimer_dg),

    "-> FORWARD PRIMER (5' -> 3'):",
    sprintf("   Sequence : %s", bp$forward_seq),
    sprintf("   Tm       : %.1f \u00B0C", bp$forward_tm),
    sprintf("   GC       : %.1f %%\n", bp$forward_gc),

    "-> REVERSE PRIMER (5' -> 3'):",
    sprintf("   Sequence : %s", bp$reverse_seq),
    sprintf("   Tm       : %.1f \u00B0C", bp$reverse_tm),
    sprintf("   GC       : %.1f %%\n", bp$reverse_gc),

    "===========================================================",
    "--- SUMMARY OF THE BEST CANDIDATES ---",
    " (Ranked from best to worst by overall score)\n"
  )

  summary_table <- utils::capture.output(
    print(design$top[, c(
      "forward_seq", "forward_tm", "reverse_seq", "reverse_tm",
      "amplicon_bp", "fragment_1", "fragment_2", "score", "recommended"
    )], row.names = FALSE, right = FALSE)
  )

  report <- c(
    report, summary_table,
    "\n===========================================================",
    "--- HOW TO READ THIS REPORT ---",
    "* Recommended: the detailed pair has the best balance of Tm, GC,",
    "  Tm difference between primers, and low (hetero)dimer risk.",
    "* Estimated RFLP fragments: computed using the SNP coordinate as an",
    "  approximation of the cut site; verify the exact value after",
    "  simulating the PCR and locating the enzyme site on the real",
    "  amplicon (see simulate_pcr() and find_restriction_site()).",
    "* Tm: calculated with the Nearest-Neighbor method (SantaLucia 2004)",
    "  and Owczarzy 2004 salt correction (see calc_tm()).",
    "* Delta G (dimer/hairpin): heuristic complementarity estimate, not a",
    "  rigorous thermodynamic calculation. Values closer to 0 are better.",
    "* Score: internal ranking index (lower value = better candidate).",
    "",
    "NOTE: before ordering synthesis, verifying the pair with a",
    "specialized tool such as IDT OligoAnalyzer is recommended."
  )

  writeLines(report, con = output_file)
  message(sprintf("Results exported to: %s", output_file))
  invisible(output_file)
}
