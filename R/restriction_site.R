#' Find a restriction enzyme site in the amplicon
#'
#' Searches for a restriction enzyme recognition motif (e.g. `"GANTC"` for
#' HinfI) within the amplicon sequence, allowing IUPAC ambiguous bases
#' (`N`, `R`, `Y`, etc.).
#'
#' @param amplicon `DNAString` object with the amplicon sequence (see
#'   [simulate_pcr()]), or the full `"rflp_amplicon"` object returned by
#'   that function.
#' @param enzyme_motif Character string with the recognition motif, in
#'   IUPAC code (default `"GANTC"`, the HinfI site).
#' @param cut_offset Position of the cut within the motif, counted from the
#'   motif's first nucleotide (1-based). Default `1`, which corresponds to
#'   HinfI's cut pattern (`G^ANTC`, cutting between positions 1 and 2 of
#'   the motif). Adjust according to the enzyme used (consult REBASE or the
#'   supplier's technical sheet).
#'
#' @return A list with:
#'   \item{sites}{`data.frame` with columns `start`, `end` and `cut_pos`
#'     (coordinates within the amplicon, 1-based).}
#'   \item{n_sites}{Number of sites found.}
#'
#' @export
find_restriction_site <- function(amplicon, enzyme_motif = "GANTC", cut_offset = 1) {
  if (inherits(amplicon, "rflp_amplicon")) amplicon <- amplicon$amplicon

  motif <- Biostrings::DNAString(enzyme_motif)
  hits <- Biostrings::matchPattern(motif, amplicon, fixed = FALSE)

  n <- length(hits)
  if (n == 0) {
    message(sprintf("No sites found for motif '%s' in the amplicon.", enzyme_motif))
    return(list(
      sites = data.frame(start = integer(), end = integer(), cut_pos = integer()),
      n_sites = 0
    ))
  }

  sites <- data.frame(start = Biostrings::start(hits), end = Biostrings::end(hits))
  sites$cut_pos <- sites$start + cut_offset - 1

  message(sprintf("Found %d site(s) for '%s' in the amplicon.", n, enzyme_motif))
  list(sites = sites, n_sites = n)
}
