#' Resolve the alternate allele for PIRA-PCR design
#'
#' Takes the reference base already identified locally by
#' \code{\link{locate_snp}} (via flanking-sequence matching against the
#' loaded FASTA) and cross-checks it against the alleles reported by dbSNP
#' (through \code{\link{fetch_dbsnp_alleles}}) for a given rsID. Returns the
#' alternate allele(s) needed to design PIRA-PCR mismatch primers.
#'
#' @param snp_base Character. Single-nucleotide base observed at the SNP
#'   position in the locally loaded sequence (forward strand), as returned
#'   by \code{locate_snp()}.
#' @param rsid Character. RefSNP identifier, e.g. \code{"rs1801133"}.
#' @param species Character. Passed to \code{\link{fetch_dbsnp_alleles}}.
#'   Default \code{"human"}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{reference_allele}{Character. The base matched locally (\code{snp_base}, as provided).}
#'   \item{alternate_alleles}{Character vector. Alternate allele(s), already
#'     expressed on the same strand as the local sequence (i.e. ready to use
#'     directly against \code{gene_seq}).}
#'   \item{strand_used}{Character. \code{"forward"} if \code{snp_base} matched
#'     an Ensembl allele directly, or \code{"complement"} if a strand flip was
#'     needed to find a match.}
#'   \item{multiallelic}{Logical. \code{TRUE} if more than one alternate
#'     allele was found (the caller must choose which one to design against).}
#'   \item{orientation_verified}{Logical. \code{TRUE} if the forward/complement
#'     strand determination was confirmed against dbSNP's ancestral allele;
#'     \code{FALSE} if it fell back to plain set-membership matching (only
#'     possible to distinguish when the locus has more than 2 alleles). See
#'     Details.}
#'   \item{dbsnp}{The full result of \code{fetch_dbsnp_alleles()}, for reference.}
#' }
#'
#' @details
#' Strand orientation is resolved by comparing \code{snp_base} against
#' dbSNP's \code{ancestral_allele} first (forward if equal, complement if
#' the complement of \code{snp_base} is equal). This ancestral-based check
#' takes priority because at multiallelic loci, \code{snp_base} can coincide
#' with a non-reference allele by chance, producing a plausible-looking but
#' incorrect direct ("forward") match. This was confirmed against real data
#' for \code{rs1801133} (MTHFR C677T), which Ensembl reports as triallelic
#' (\code{G/A/C}) because of two distinct clinically documented substitutions
#' at the same nucleotide (c.665C>T and c.665C>G). Naive set-membership
#' matching alone selects the wrong alternate allele set for that locus;
#' the ancestral-based check resolves it correctly. If no ancestral allele
#' is available, or \code{snp_base} is a valid non-ancestral reference base,
#' the function falls back to plain set-membership matching and flags the
#' result via \code{orientation_verified = FALSE} when the locus is
#' multiallelic.
#'
#' If \code{snp_base} cannot be matched to any dbSNP allele (forward or
#' complement), the function issues a \code{warning} and returns
#' \code{alternate_alleles = character(0)} with \code{strand_used = NA} —
#' it does not stop execution, since this may reflect a genuinely
#' informative discrepancy (assembly mismatch, mislocalized SNP) that the
#' user should investigate rather than have silently masked.
#'
#' @examples
#' \dontrun{
#' # MTHFR C677T: reference C, alternate T
#' resolve_pira_alleles(snp_base = "C", rsid = "rs1801133")
#' }
#'
#' @export
resolve_pira_alleles <- function(snp_base, rsid, species = "human") {

  if (missing(snp_base) || !is.character(snp_base) || length(snp_base) != 1 ||
      !nzchar(snp_base)) {
    stop("'snp_base' must be a single-character string, e.g. \"C\".", call. = FALSE)
  }
  snp_base <- toupper(snp_base)
  if (!snp_base %in% c("A", "C", "G", "T")) {
    stop("'snp_base' must be one of A, C, G, T. Got: '", snp_base, "'.",
         call. = FALSE)
  }

  dbsnp <- fetch_dbsnp_alleles(rsid = rsid, species = species)

  complement_base <- function(b) {
    map <- c(A = "T", T = "A", C = "G", G = "C")
    unname(map[b])
  }

  alleles <- dbsnp$alleles
  ancestral <- dbsnp$ancestral_allele

  # --- Determine strand orientation -----------------------------------------
  # Preferred method: compare snp_base against the dbSNP ancestral allele.
  # This resolves cases where the locus has more than 2 observed alleles and
  # snp_base happens to coincide with a *non-reference* allele on the wrong
  # strand interpretation purely by chance (verified against real data for
  # rs1801133 / MTHFR C677T, which is genuinely triallelic: G/A/C on the
  # Ensembl + strand, corresponding to two distinct clinical substitutions,
  # c.665C>T and c.665C>G, on the coding strand). Naive "is snp_base anywhere
  # in the allele set" matching cannot distinguish these cases reliably.
  orientation <- NA_character_
  ancestral_available <- !is.na(ancestral) && ancestral %in% c("A", "C", "G", "T")

  if (ancestral_available) {
    if (identical(snp_base, ancestral)) {
      orientation <- "forward"
    } else if (identical(complement_base(snp_base), ancestral)) {
      orientation <- "complement"
    }
  }

  ancestral_based <- !is.na(orientation)

  if (is.na(orientation)) {
    # Fallback: no ancestral allele available, or snp_base is a valid
    # non-ancestral (derived) reference base that doesn't match the
    # ancestral allele on either strand. Fall back to plain set membership,
    # as before, but flag the result as strand-unverified when the locus
    # has more than 2 alleles (the ambiguous scenario this fallback cannot
    # fully resolve).
    if (snp_base %in% alleles) {
      orientation <- "forward"
    } else if (complement_base(snp_base) %in% alleles) {
      orientation <- "complement"
    }
  }

  if (is.na(orientation)) {
    warning(
      "Local reference base ('", snp_base, "') does not match any dbSNP ",
      "allele (", paste(alleles, collapse = "/"), ") for '", rsid, "', even ",
      "after checking the complementary strand. This may indicate an ",
      "assembly mismatch, an incorrect rsID, or a mislocalized SNP position. ",
      "Returning no alternate allele; please investigate before proceeding.",
      call. = FALSE
    )
    return(list(
      reference_allele  = snp_base,
      alternate_alleles = character(0),
      strand_used       = NA_character_,
      multiallelic       = NA,
      orientation_verified = FALSE,
      dbsnp             = dbsnp
    ))
  }

  strand_used <- orientation
  matched_allele <- if (orientation == "forward") snp_base else complement_base(snp_base)
  alt_alleles_dbsnp <- setdiff(alleles, matched_allele)
  alt_alleles_local <- if (orientation == "forward") {
    alt_alleles_dbsnp
  } else {
    vapply(alt_alleles_dbsnp, complement_base, character(1))
  }

  if (orientation == "complement") {
    warning(
      "Local reference base ('", snp_base, "') did not match dbSNP alleles ",
      "directly (", paste(alleles, collapse = "/"), ") for '", rsid, "', but ",
      "its complement did", if (ancestral_based) " (confirmed against the ancestral allele)" else "",
      ". Assuming the local sequence is on the opposite strand from the ",
      "Ensembl mapping; alternate allele(s) have been converted back to the ",
      "local strand. Please verify this assumption.",
      call. = FALSE
    )
  }

  multiallelic <- length(alt_alleles_local) > 1
  if (multiallelic) {
    warning(
      "rsID '", rsid, "' is multiallelic: ", length(alt_alleles_local),
      " alternate allele(s) found (", paste(alt_alleles_local, collapse = ", "),
      "). PIRA-PCR primer design must target one specific alternate allele; ",
      "this function does not choose automatically.", call. = FALSE
    )
  }
  if (multiallelic && !ancestral_based) {
    warning(
      "Strand orientation for '", rsid, "' could not be confirmed against ",
      "the ancestral allele (unavailable, or snp_base is a non-ancestral ",
      "reference base). With more than 2 alleles present, the alternate ",
      "allele(s) reported here should be treated as provisional until ",
      "manually verified.", call. = FALSE
    )
  }

  list(
    reference_allele  = snp_base,
    alternate_alleles = alt_alleles_local,
    strand_used       = strand_used,
    multiallelic       = multiallelic,
    orientation_verified = ancestral_based,
    dbsnp             = dbsnp
  )
}
