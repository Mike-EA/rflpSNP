#' Run the complete tetra-primer ARMS-PCR workflow
#'
#' @param gene_seq Reference `Biostrings::DNAString`.
#' @param flank_seq SNP flanking sequence for [locate_snp()].
#' @param alt_allele Explicit alternative base.
#' @param export_txt Whether to write the text report.
#' @param output_file Destination for the optional report.
#' @param ... Additional arguments for [design_arms_primers()].
#'
#' @return A list containing the SNP location, design, optional report path,
#'   simulated products and gel plot.
#' @export
run_arms_pcr_pipeline <- function(gene_seq, flank_seq, alt_allele, export_txt = TRUE,
                                  output_file = "arms_pcr_primers_results.txt", ...) {
  snp <- locate_snp(gene_seq, flank_seq)
  design <- design_arms_primers(gene_seq, snp$snp_pos, snp$snp_base, alt_allele, ...)
  if (!design$n_valid_sets) stop("No valid ARMS-PCR set was found; inspect design$diagnostics and design$exclusion_diagnostics.")
  report <- if (export_txt) export_arms_primers_txt(design, output_file) else NULL
  pcr <- simulate_arms_pcr(gene_seq, design$best_set, snp$snp_pos, snp$snp_base, alt_allele)
  list(snp = snp, design = design, report = report, pcr_result = pcr, gel = simulate_arms_gel(pcr))
}
