#' Locate a SNP using its flanking sequence
#'
#' Locates the coordinate of a SNP within a reference gene sequence using
#' the 5' flanking sequence reported by dbSNP (NIH) as an anchor. It first
#' searches the forward strand; if there is no match, it searches the
#' complementary strand and translates the coordinate back to the forward
#' reference frame.
#'
#' @param gene_seq `Biostrings::DNAString` object with the complete gene
#'   sequence (see [read_gene_fasta()]).
#' @param flank_seq Character string or `DNAString` with the 5' flanking
#'   sequence obtained from dbSNP, immediately preceding the SNP position.
#' @param snp_offset Integer. Number of bases after the end of the flanking
#'   match where the SNP is located (default `1`, i.e. the base immediately
#'   following the flank). Adjust only if the flanking-sequence convention
#'   you are using differs from dbSNP's.
#' @param rsid Character. Optional RefSNP identifier (e.g. `"rs1801133"`).
#'   When provided, `locate_snp()` additionally calls
#'   [resolve_pira_alleles()] to fetch and cross-check the alternate
#'   allele(s) needed for PIRA-PCR primer design (requires internet access).
#'   When `NULL` (the default), behavior is unchanged from previous versions
#'   of this function.
#' @param species Character. Passed to [resolve_pira_alleles()] when `rsid`
#'   is provided. Default `"human"`.
#'
#' @return A list with the elements:
#'   \item{snp_pos}{1-based coordinate of the SNP in `gene_seq`.}
#'   \item{snp_base}{Nucleotide found at `snp_pos`.}
#'   \item{strand}{`"forward"` or `"complementary"`, the strand on which the
#'     flanking sequence match was found.}
#'   \item{alleles}{Only present when `rsid` is supplied: the full result of
#'     [resolve_pira_alleles()] (reference/alternate allele(s), strand
#'     orientation, multiallelic flag, etc.).}
#'
#' @details
#' This function assumes there is exactly one relevant match of the
#' flanking sequence in the reference sequence. If the flanking sequence
#' appears more than once, the first match is used and a warning is
#' issued; in that case, using a longer and more specific flanking sequence
#' is recommended.
#'
#' When `rsid` is supplied, allele resolution errors (e.g. no internet
#' connection, invalid rsID) propagate as errors from `locate_snp()` itself
#' — the SNP has already been located successfully by that point, so this
#' distinguishes "the SNP location logic failed" from "the SNP was located,
#' but dbSNP allele resolution failed," while still stopping the pipeline
#' either way (a PIRA-PCR design cannot proceed without a resolved
#' alternate allele).
#'
#' @examples
#' \dontrun{
#' result <- locate_snp(gene_seq, flank_seq, rsid = "rs1801133")
#' result$alleles$alternate_alleles
#' }
#'
#' @export
locate_snp <- function(gene_seq, flank_seq, snp_offset = 1,
                        rsid = NULL, species = "human") {
  if (is.character(flank_seq)) {
    flank_seq <- Biostrings::DNAString(flank_seq)
  }

  fwd_matches <- Biostrings::matchPattern(flank_seq, gene_seq, fixed = TRUE)

  if (length(fwd_matches) > 0) {
    if (length(fwd_matches) > 1) {
      warning(sprintf(
        "The flanking sequence matches %d times on the forward strand; using the first match.",
        length(fwd_matches)
      ))
    }
    fwd_end <- Biostrings::end(fwd_matches)[1]
    snp_pos <- fwd_end + snp_offset

    if (snp_pos > length(gene_seq) || snp_pos < 1) {
      stop("The computed SNP position falls outside the sequence bounds (forward strand). Check snp_offset and the flanking sequence.")
    }

    snp_base <- as.character(gene_seq[snp_pos])
    strand_info <- "forward"
  } else {
    rev_seq <- Biostrings::reverseComplement(gene_seq)
    rev_matches <- Biostrings::matchPattern(flank_seq, rev_seq, fixed = TRUE)

    if (length(rev_matches) == 0) {
      stop("The flanking sequence was not found on either strand. Verify that it belongs to the loaded gene.")
    }
    if (length(rev_matches) > 1) {
      warning(sprintf(
        "The flanking sequence matches %d times on the complementary strand; using the first match.",
        length(rev_matches)
      ))
    }

    rev_end <- Biostrings::end(rev_matches)[1]
    snp_pos <- length(gene_seq) - rev_end - (snp_offset - 1)

    if (snp_pos < 1 || snp_pos > length(gene_seq)) {
      stop("The computed SNP position falls outside the sequence bounds (complementary strand). Check snp_offset and the flanking sequence.")
    }

    snp_base <- as.character(gene_seq[snp_pos])
    strand_info <- "complementary"
  }

  message(sprintf(
    "Nucleotide at position %d (%s) [found on the %s strand]",
    snp_pos, snp_base, strand_info
  ))

  result <- list(snp_pos = snp_pos, snp_base = snp_base, strand = strand_info)

  if (!is.null(rsid)) {
    result$alleles <- resolve_pira_alleles(
      snp_base = snp_base, rsid = rsid, species = species
    )
  }

  result
}
