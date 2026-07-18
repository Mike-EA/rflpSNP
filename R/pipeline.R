#' Run the complete end-to-end PCR-RFLP design workflow
#'
#' Convenience function that chains [locate_snp()], [design_primers()],
#' [simulate_pcr()], [find_restriction_site()] and, optionally,
#' [export_primers_txt()], using the best primer pair found. Intended for a
#' quick first pass; for teaching purposes, also running each function
#' separately is recommended so that students understand each stage of the
#' design.
#'
#' @param gene_seq `DNAString` object with the complete gene sequence (see
#'   [read_gene_fasta()]).
#' @param flank_seq dbSNP flanking sequence (see [locate_snp()]).
#' @param enzyme_motif Restriction enzyme recognition motif (see
#'   [find_restriction_site()]).
#' @param cut_offset Cut position within the motif (see
#'   [find_restriction_site()]).
#' @param export_txt If `TRUE` (default), exports the primer report to text
#'   via [export_primers_txt()].
#' @param output_file Path of the text file to export (if
#'   `export_txt = TRUE`).
#' @param ... Additional arguments passed to [design_primers()] (e.g.
#'   `tm_min`, `amplicon_max`, `min_fragment_diff`, etc.).
#'
#' @return A list with the elements `snp`, `design`, `pcr_result` and
#'   `restriction_result`, i.e. the result of each intermediate stage.
#'
#' @export
run_pcr_rflp_pipeline <- function(gene_seq, flank_seq,
                                   enzyme_motif = "GANTC", cut_offset = 1,
                                   export_txt = TRUE,
                                   output_file = "primers_PCR_RFLP_results.txt",
                                   ...) {
  snp <- locate_snp(gene_seq, flank_seq)

  design <- design_primers(gene_seq, snp_pos = snp$snp_pos, ...)

  if (export_txt) {
    export_primers_txt(design, output_file = output_file)
  }

  bp <- design$best_pair
  pcr_result <- simulate_pcr(gene_seq, bp$forward_seq, bp$reverse_seq)
  restriction_result <- find_restriction_site(pcr_result, enzyme_motif = enzyme_motif, cut_offset = cut_offset)

  if (restriction_result$n_sites == 0) {
    warning(
      "The best primer pair does not produce an amplicon containing the ",
      "given restriction site. This can happen if the enzyme motif does ",
      "not match the reference allele in the loaded FASTA (i.e. the site ",
      "only exists on the alternate allele). Consider reviewing other ",
      "pairs in design$top, or generating the alternate-allele amplicon ",
      "by substituting the SNP base before calling simulate_pcr()."
    )
  }

  list(
    snp = snp,
    design = design,
    pcr_result = pcr_result,
    restriction_result = restriction_result
  )
}
