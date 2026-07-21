#' IUPAC nucleotide ambiguity code table
#'
#' Internal lookup table mapping each IUPAC nucleotide code to the literal
#' base(s) it represents.
#' @noRd
.IUPAC_TABLE <- list(
  A = "A", C = "C", G = "G", T = "T",
  R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"),
  B = c("C", "G", "T"), D = c("A", "G", "T"), H = c("A", "C", "T"), V = c("A", "C", "G"),
  N = c("A", "C", "G", "T")
)

#' Literal bases represented by an IUPAC code
#'
#' @param code Single-character IUPAC nucleotide code.
#' @return Character vector of literal bases (A/C/G/T) matching `code`.
#' @export
iupac_bases <- function(code) {
  code <- toupper(code)
  bases <- .IUPAC_TABLE[[code]]
  if (is.null(bases)) {
    stop("Unrecognized IUPAC code: '", code, "'.", call. = FALSE)
  }
  bases
}

#' Test whether a literal base satisfies an IUPAC code
#'
#' @param base Single literal nucleotide (A/C/G/T).
#' @param code Single IUPAC nucleotide code (may itself be A/C/G/T, or an
#'   ambiguity code such as R, Y, N, etc.).
#' @return Logical.
#' @export
iupac_matches <- function(base, code) {
  base <- toupper(base)
  toupper(base) %in% iupac_bases(code)
}
