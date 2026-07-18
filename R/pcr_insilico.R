#' Simulate an in silico PCR
#'
#' Given a reference sequence and a primer pair, locates the binding sites
#' (allowing a configurable number of mismatches to assess specificity) and
#' extracts the resulting amplicon as a `DNAString` object.
#'
#' @param gene_seq `DNAString` object, the complete reference sequence (not
#'   necessarily the same variable used in [design_primers()]; it can be
#'   any sequence on which you want to test the primer pair).
#' @param fwd_primer Character string or `DNAString` with the forward
#'   primer (5' -> 3').
#' @param rev_primer Character string or `DNAString` with the reverse
#'   primer (5' -> 3').
#' @param max_mismatch Maximum number of mismatches allowed when searching
#'   for additional binding sites, to assess specificity (default `2`). The
#'   amplicon itself is always computed from the exact match.
#'
#' @return A list of class `"rflp_amplicon"` with:
#'   \item{amplicon}{`DNAString` object with the amplicon sequence.}
#'   \item{start, end}{1-based coordinates in `gene_seq`.}
#'   \item{size}{Amplicon size in bp.}
#'   \item{n_sites_fwd, n_sites_rev}{Number of binding sites found allowing
#'     `max_mismatch` mismatches; values greater than 1 suggest possible
#'     non-specific binding in the reference sequence used.}
#'
#' @export
simulate_pcr <- function(gene_seq, fwd_primer, rev_primer, max_mismatch = 2) {
  if (is.character(fwd_primer)) fwd_primer <- Biostrings::DNAString(fwd_primer)
  if (is.character(rev_primer)) rev_primer <- Biostrings::DNAString(rev_primer)

  rev_primer_rc <- Biostrings::reverseComplement(rev_primer)

  match_fwd_exact <- Biostrings::matchPattern(fwd_primer, gene_seq)
  match_rev_exact <- Biostrings::matchPattern(rev_primer_rc, gene_seq)

  if (length(match_fwd_exact) == 0 || length(match_rev_exact) == 0) {
    stop("No exact match was found for one or both primers in the provided reference sequence.")
  }

  start_amp <- Biostrings::start(match_fwd_exact)[1]
  end_amp <- Biostrings::end(match_rev_exact)[1]

  if (end_amp <= start_amp) {
    stop("The reverse primer position precedes or coincides with the forward primer's; check the orientation and sequences of the primers.")
  }

  size_amp <- end_amp - start_amp + 1
  amplicon <- Biostrings::subseq(gene_seq, start_amp, end_amp)

  match_fwd_mm <- Biostrings::matchPattern(fwd_primer, gene_seq, max.mismatch = max_mismatch, fixed = FALSE)
  match_rev_mm <- Biostrings::matchPattern(rev_primer_rc, gene_seq, max.mismatch = max_mismatch, fixed = FALSE)

  message(sprintf("Amplicon: %d bp (positions %d-%d)", size_amp, start_amp, end_amp))
  message(sprintf(
    "Binding sites (<=%d mismatches) - Fwd: %d | Rev: %d",
    max_mismatch, length(match_fwd_mm), length(match_rev_mm)
  ))
  if (length(match_fwd_mm) > 1 || length(match_rev_mm) > 1) {
    warning("Multiple potential binding sites were detected; review primer specificity before continuing.")
  }

  result <- list(
    amplicon = amplicon, start = start_amp, end = end_amp, size = size_amp,
    n_sites_fwd = length(match_fwd_mm), n_sites_rev = length(match_rev_mm)
  )
  class(result) <- "rflp_amplicon"
  result
}

#' @export
print.rflp_amplicon <- function(x, ...) {
  cat("<rflp_amplicon>\n")
  cat(sprintf("Size: %d bp (positions %d-%d in the reference sequence)\n", x$size, x$start, x$end))
  cat(sprintf("Potential binding sites - Fwd: %d | Rev: %d\n", x$n_sites_fwd, x$n_sites_rev))
  cat("Sequence:\n")
  print(x$amplicon)
  invisible(x)
}
