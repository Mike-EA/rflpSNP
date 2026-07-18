#' Define the working region around the SNP
#'
#' Given the SNP coordinate, delimits the subregion of the gene sequence
#' within which candidate primers will be searched. The region is built
#' **asymmetrically** with respect to the SNP (different lengths upstream
#' and downstream) so that, after digestion with the restriction enzyme,
#' the two resulting fragments have clearly distinguishable sizes on an
#' agarose gel (the whole point of a PCR-RFLP assay is being able to tell
#' genotypes apart by their band pattern).
#'
#' @param gene_seq `DNAString` object with the complete gene sequence.
#' @param snp_pos 1-based coordinate of the SNP within `gene_seq`, as
#'   returned by [locate_snp()].
#' @param upstream Number of bp to include upstream of the SNP (default
#'   `160`).
#' @param downstream Number of bp to include downstream of the SNP (default
#'   `200`).
#' @param max_region_size Recommended maximum size (bp) for the working
#'   region in a standard PCR-RFLP assay (default `500`). A warning is
#'   issued if exceeded; execution is not stopped, since in some cases
#'   (genes with few primer options) it may be necessary.
#'
#' @return A list with:
#'   \item{start_region, end_region}{1-based coordinates in `gene_seq`.}
#'   \item{sequence_region}{`DNAString` object of the subregion.}
#'   \item{snp_pos}{The SNP coordinate, for convenience in downstream
#'     functions.}
#'
#' @export
define_work_region <- function(gene_seq, snp_pos,
                                upstream = 160, downstream = 200,
                                max_region_size = 500) {
  if (upstream == downstream) {
    warning(
      "upstream and downstream are equal: the SNP will end up centered in ",
      "the working region, which makes it harder to obtain RFLP fragments ",
      "of distinguishable sizes on a gel. Using different values is ",
      "recommended (e.g. upstream = 160, downstream = 200)."
    )
  }

  region_size <- upstream + downstream
  if (region_size > max_region_size) {
    warning(sprintf(
      "The requested working region (%d bp) exceeds the recommended maximum for PCR-RFLP (%d bp).",
      region_size, max_region_size
    ))
  }

  start_region <- max(1, snp_pos - upstream)
  end_region <- min(length(gene_seq), snp_pos + downstream)

  if (end_region <= start_region) {
    stop("The resulting working region is invalid (end_region <= start_region). Check snp_pos, upstream and downstream.")
  }

  sequence_region <- Biostrings::subseq(gene_seq, start_region, end_region)

  message(sprintf(
    "Working region: nucleotide %d to %d (length: %d bp; SNP at relative position %d)",
    start_region, end_region, length(sequence_region), snp_pos - start_region + 1
  ))

  list(
    start_region = start_region,
    end_region = end_region,
    sequence_region = sequence_region,
    snp_pos = snp_pos
  )
}
