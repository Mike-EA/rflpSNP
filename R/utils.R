#' rflpSNP: In Silico Design and Simulation of SNP Genotyping Assays
#'
#' The \pkg{rflpSNP} package automates the complete workflow for designing a
#' PCR-RFLP or tetra-primer ARMS-PCR assay for genotyping a known SNP:
#'
#' \enumerate{
#'   \item Load a reference sequence (\code{\link{read_gene_fasta}}).
#'   \item Locate the SNP using its dbSNP flanking sequence
#'     (\code{\link{locate_snp}}).
#'   \item Define the working region around the SNP
#'     (\code{\link{define_work_region}}).
#'   \item Design and filter candidate primer pairs
#'     (\code{\link{design_primers}}).
#'   \item Export the results to text (\code{\link{export_primers_txt}}).
#'   \item Simulate the PCR in silico (\code{\link{simulate_pcr}}).
#'   \item Locate the restriction enzyme cut site
#'     (\code{\link{find_restriction_site}}).
#'   \item Plot the amplicon and the highlighted sequence
#'     (\code{\link{plot_amplicon_map}}, \code{\link{plot_sequence_map}}).
#'   \item Simulate the expected agarose gel (\code{\link{simulate_gel}}).
#' }
#'
#' The function \code{\link{run_pcr_rflp_pipeline}} chains all of these steps
#' together for a quick end-to-end workflow.
#' For SNPs without a useful restriction-site change,
#' \code{\link{run_arms_pcr_pipeline}} designs four ARMS primers and returns
#' the report, genotype products and virtual gel.
#'
#' @keywords internal
"_PACKAGE"

# Variables used inside ggplot2::aes()/aggregate() that are not global R
# variables but column names evaluated through non-standard evaluation
# (NSE). They are declared here purely so that R CMD check does not report
# false positives of the "no visible binding for global variable" kind.
utils::globalVariables(c(
  "Column", "Row", "Type", "Base", "Label", "GenPos",
  "Size", "Condition", "Y_pos", "Intensity", "Lane", "X",
  "lane", "y", "product", "size_bp"
))
