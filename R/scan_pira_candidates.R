#' Scan for PIRA-PCR restriction site candidates
#'
#' For a SNP already localized in `gene_seq` (see [locate_snp()]), searches
#' the restriction enzyme catalog (see [restriction_enzymes]) for enzymes
#' whose recognition site can be completed for exactly one allele (reference
#' or alternate) by introducing 1 or 2 artificial mismatches in a primer,
#' within the window defined by [PIRA_MISMATCH_WINDOW] relative to the SNP
#' (which sits at the primer's 3' terminus, ARMS-style).
#'
#' @param gene_seq `Biostrings::DNAString` with the reference sequence.
#' @param snp_pos Integer. 1-based position of the SNP in `gene_seq` (see
#'   [locate_snp()]).
#' @param reference_allele,alternate_allele Single-character strings (A/C/G/T),
#'   in `gene_seq`'s own forward-strand coordinate system (see
#'   [resolve_pira_alleles()]).
#' @param enzymes A data frame with the same structure as
#'   [restriction_enzymes] (the default). Must include at least `enzyme`,
#'   `recognition_site`, and `cut_type`.
#' @param mismatch_window Integer vector of allowed artificial-mismatch
#'   positions, counted from the primer's 3' end (the SNP position itself is
#'   position 0 and is not part of this window). Default
#'   [PIRA_MISMATCH_WINDOW] (`1:3`).
#' @param max_artificial_mismatches Integer. Maximum number of artificial
#'   mismatches allowed in a single design. Default `2`.
#' @param orientations Character vector, any of `"forward"`, `"reverse"`.
#'   Default both.
#'
#' @return A data frame, one row per valid candidate, with columns:
#' \describe{
#'   \item{enzyme}{Enzyme name.}
#'   \item{orientation}{`"forward"` or `"reverse"`.}
#'   \item{site_start, site_end}{Recognition site coordinates in `gene_seq`.}
#'   \item{n_artificial_mismatches}{1 or 2.}
#'   \item{mismatch_positions}{Semicolon-separated genomic positions of the
#'     artificial mismatches.}
#'   \item{mismatch_bases_required}{Semicolon-separated list of the base(s)
#'     that satisfy the enzyme's pattern at each mismatch position (a
#'     position may have more than one valid base when the pattern is
#'     itself an IUPAC ambiguity code; listed as e.g. `"A/G"`).}
#'   \item{creates_site_for}{`"reference"` or `"alternate"` -- the allele
#'     for which the site is completed.}
#'   \item{cut_type, recognition_site}{Carried over from the enzyme catalog.}
#' }
#' Returns a zero-row data frame (same columns) if no candidates are found.
#'
#' @details
#' The comparison is allele-vs-allele on the same template and the same
#' choice of artificial-mismatch bases (see the Fase 3 planning notes in the
#' package source for the full derivation). A design is only reported if
#' the site is completed for exactly one of the two alleles -- if both or
#' neither allele completes it, the site would not discriminate genotypes
#' and is not reported.
#'
#' Coordinates are always in `gene_seq`'s forward-strand frame, for both
#' orientations. For `"reverse"`, the mismatch window falls downstream of
#' `snp_pos` (higher coordinates); constructing the actual reverse-primer
#' oligo sequence (reverse-complementing this region) is left to a later
#' stage of the pipeline.
#'
#' @examples
#' \dontrun{
#' scan_pira_candidates(gene_seq, snp_pos = 120,
#'                       reference_allele = "C", alternate_allele = "A")
#' }
#'
#' @export
scan_pira_candidates <- function(gene_seq, snp_pos,
                                  reference_allele, alternate_allele,
                                  enzymes = restriction_enzymes,
                                  mismatch_window = PIRA_MISMATCH_WINDOW,
                                  max_artificial_mismatches = 2,
                                  orientations = c("forward", "reverse")) {

  # --- Validation ------------------------------------------------------
  if (!methods::is(gene_seq, "DNAString")) {
    stop("'gene_seq' must be a Biostrings::DNAString object.", call. = FALSE)
  }
  seq_len <- length(gene_seq)
  if (!is.numeric(snp_pos) || length(snp_pos) != 1 || snp_pos < 1 || snp_pos > seq_len) {
    stop("'snp_pos' must be a single integer within the bounds of 'gene_seq'.", call. = FALSE)
  }
  snp_pos <- as.integer(snp_pos)

  reference_allele <- toupper(reference_allele)
  alternate_allele <- toupper(alternate_allele)
  if (!reference_allele %in% c("A", "C", "G", "T") ||
      !alternate_allele %in% c("A", "C", "G", "T")) {
    stop("'reference_allele' and 'alternate_allele' must each be one of A, C, G, T.", call. = FALSE)
  }
  if (identical(reference_allele, alternate_allele)) {
    stop("'reference_allele' and 'alternate_allele' must differ.", call. = FALSE)
  }

  orientations <- match.arg(orientations, c("forward", "reverse"), several.ok = TRUE)

  required_cols <- c("enzyme", "recognition_site", "cut_type")
  missing_cols <- setdiff(required_cols, names(enzymes))
  if (length(missing_cols) > 0) {
    stop("'enzymes' is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  empty_result <- data.frame(
    enzyme = character(0), orientation = character(0),
    site_start = integer(0), site_end = integer(0),
    n_artificial_mismatches = integer(0),
    mismatch_positions = character(0), mismatch_bases_required = character(0),
    creates_site_for = character(0),
    cut_type = character(0), recognition_site = character(0),
    stringsAsFactors = FALSE
  )

  results <- vector("list", 0)

  for (i in seq_len(nrow(enzymes))) {
    site <- toupper(enzymes$recognition_site[i])
    L <- nchar(site)
    if (L < 1 || L > seq_len) next
    site_chars <- strsplit(site, "")[[1]]

    for (orientation in orientations) {
      window_positions <- if (orientation == "forward") {
        snp_pos - mismatch_window
      } else {
        snp_pos + mismatch_window
      }
      window_positions <- window_positions[window_positions >= 1 & window_positions <= seq_len]
      if (length(window_positions) == 0) next

      site_start_range <- (snp_pos - L + 1):snp_pos

      for (site_start in site_start_range) {
        site_end <- site_start + L - 1
        if (site_start < 1 || site_end > seq_len) next
        if (!any(window_positions >= site_start & window_positions <= site_end)) next

        ok <- TRUE
        n_mismatch <- 0L
        mismatch_positions <- integer(0)
        mismatch_bases <- character(0)
        code_at_snp <- NA_character_

        for (offset in seq_len(L)) {
          pos <- site_start + offset - 1L
          code <- site_chars[offset]

          if (pos == snp_pos) {
            code_at_snp <- code
            next
          }

          genomic_base <- as.character(gene_seq[pos])

          if (pos %in% window_positions) {
            if (!iupac_matches(genomic_base, code)) {
              n_mismatch <- n_mismatch + 1L
              mismatch_positions <- c(mismatch_positions, pos)
              mismatch_bases <- c(mismatch_bases, paste(iupac_bases(code), collapse = "/"))
            }
          } else {
            if (!iupac_matches(genomic_base, code)) {
              ok <- FALSE
              break
            }
          }
        }

        if (!ok) next
        if (is.na(code_at_snp)) next  # defensive: site didn't actually cover snp_pos
        if (n_mismatch < 1 || n_mismatch > max_artificial_mismatches) next

        ref_match <- iupac_matches(reference_allele, code_at_snp)
        alt_match <- iupac_matches(alternate_allele, code_at_snp)

        if (xor(ref_match, alt_match)) {
          creates_for <- if (ref_match) "reference" else "alternate"
          results[[length(results) + 1]] <- data.frame(
            enzyme = enzymes$enzyme[i],
            orientation = orientation,
            site_start = site_start,
            site_end = site_end,
            n_artificial_mismatches = n_mismatch,
            mismatch_positions = paste(mismatch_positions, collapse = ";"),
            mismatch_bases_required = paste(mismatch_bases, collapse = ";"),
            creates_site_for = creates_for,
            cut_type = enzymes$cut_type[i],
            recognition_site = site,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(results) == 0) return(empty_result)
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
