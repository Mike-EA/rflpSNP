# Core PIRA-PCR compatibility scanner for a single (enzyme, frame) pair,
# expressed entirely in the "forward/upstream-primer" reference frame:
# the SNP always sits immediately after (or `gap` bases after) the primer
# whose 3' end is being engineered, reading 5' -> 3'. The downstream-primer
# (reverse strand) case is handled by find_pira_sites() by reverse-
# complementing the relevant flanks before calling this function, which
# keeps the indexing logic in one place.
#
# Arguments:
#   upstream_flank   character vector, 5'->3', genomic bases immediately
#                     upstream of the SNP (last element is the base right
#                     before the gap/SNP).
#   downstream_flank  character vector, 5'->3', genomic bases immediately
#                     downstream of the SNP.
#   ref_base, alt_base single-character alleles at the SNP position, on
#                     the same strand as upstream_flank/downstream_flank.
#   motif             character vector (single bases), the IUPAC motif,
#                     5'->3', length L.
#   snp_motif_pos     integer 1..L: which position of the motif aligns
#                     with the SNP.
#   min_dist_3prime   minimum allowed distance (bp) between the primer's
#                     3' terminal base and the nearest engineered mismatch
#                     (distance 1 = the terminal base itself).
#   gap               number of unmodifiable genomic bases between the
#                     primer's true 3' end and the SNP (0 = SNP is the
#                     base immediately following the primer).
#
# Returns NULL if this (enzyme, frame) combination is not viable, or a
# list describing the required engineering otherwise.
#' @keywords internal
.pira_scan_frame <- function(upstream_flank, downstream_flank, ref_base, alt_base,
                              motif, snp_motif_pos,
                              min_dist_3prime = 2, gap = 0) {
  L <- length(motif)
  s <- snp_motif_pos
  if (s < 1 || s > L) return(NULL)

  n_before <- s - 1
  n_after  <- L - s

  if (n_before > length(upstream_flank)) return(NULL)
  if (n_after > length(downstream_flank)) return(NULL)

  # --- bases downstream of the SNP: unmodifiable, must already match ---
  if (n_after > 0) {
    downstream_needed <- downstream_flank[seq_len(n_after)]
    downstream_codes <- motif[(s + 1):L]
    ok_after <- all(mapply(.iupac_matches, downstream_needed, downstream_codes))
    if (!ok_after) return(NULL)
  }

  # --- the SNP position itself must discriminate ref vs alt ---
  ref_ok <- .iupac_matches(ref_base, motif[s])
  alt_ok <- .iupac_matches(alt_base, motif[s])
  if (ref_ok == alt_ok) return(NULL)
  cuts_on <- if (ref_ok) "ref" else "alt"

  # --- bases upstream of the SNP: split into a fixed "buffer" zone
  #     (nearest the SNP, size = gap, unmodifiable) and an "engineerable"
  #     zone (farther from the SNP, where mismatches may be introduced) ---
  n_buffer <- min(gap, n_before)
  n_engineerable <- n_before - n_buffer

  if (n_buffer > 0) {
    buffer_codes <- motif[(n_engineerable + 1):n_before]
    buffer_bases <- utils::tail(upstream_flank, n_buffer)
    ok_buffer <- all(mapply(.iupac_matches, buffer_bases, buffer_codes))
    if (!ok_buffer) return(NULL)
  }

  if (n_engineerable == 0) {
    return(list(
      n_mismatch = 0, natural_site = TRUE, cuts_on = cuts_on,
      mismatch_dist_3prime = integer(0),
      required_base = character(0), original_base = character(0)
    ))
  }

  engineer_codes <- motif[seq_len(n_engineerable)]
  total_needed <- n_engineerable + n_buffer
  if (total_needed > length(upstream_flank)) return(NULL)
  relevant <- utils::tail(upstream_flank, total_needed)
  engineer_bases_actual <- relevant[seq_len(n_engineerable)]

  is_mismatch <- !mapply(.iupac_matches, engineer_bases_actual, engineer_codes)
  n_mismatch <- sum(is_mismatch)

  # Distance from the primer's 3' terminal base. Measured strictly within
  # the primer itself: the engineerable position closest to the SNP (last
  # in the vector) IS the primer's own terminal base (distance 1),
  # regardless of gap -- gap only describes unmodifiable buffer bases that
  # sit OUTSIDE the primer, between its 3' end and the SNP, and must not
  # shift this distance.
  dist_from_3prime <- rev(seq_len(n_engineerable))

  if (n_mismatch > 0) {
    min_dist_mismatch <- min(dist_from_3prime[is_mismatch])
    if (min_dist_mismatch < min_dist_3prime) return(NULL)
  }

  list(
    n_mismatch = n_mismatch,
    natural_site = (n_mismatch == 0),
    cuts_on = cuts_on,
    mismatch_dist_3prime = if (n_mismatch > 0) dist_from_3prime[is_mismatch] else integer(0),
    required_base = if (n_mismatch > 0) {
      vapply(engineer_codes[is_mismatch], function(cd) .iupac_expand(cd)[1], character(1))
    } else character(0),
    original_base = if (n_mismatch > 0) engineer_bases_actual[is_mismatch] else character(0)
  )
}

# Reverse-complements a character vector of single bases (5'->3' in,
# 5'->3' out), for mirroring the downstream/reverse-primer scenario into
# the same reference frame used by .pira_scan_frame().
#' @keywords internal
.revcomp_vec <- function(bases) {
  rev(vapply(bases, .complement_base, character(1), USE.NAMES = FALSE))
}

#' Scan for viable PIRA-PCR (enzyme, mismatch) combinations at a SNP
#'
#' For every enzyme in `enzyme_panel` and every possible alignment of its
#' recognition motif over the primer/SNP boundary, evaluates whether a
#' restriction site that discriminates the reference and alternate alleles
#' can be formed -- either naturally (0 mismatches, meaning ordinary
#' PCR-RFLP would already work and PIRA-PCR is not actually needed) or by
#' introducing a single engineered mismatch near the 3' end of one primer.
#'
#' @param gene_seq `DNAString` (or plain character) with the reference
#'   sequence.
#' @param snp_pos 1-based coordinate of the SNP in `gene_seq`.
#' @param ref_allele,alt_allele The two alleles at the SNP position
#'   (single-character strings, on the same strand as `gene_seq`), e.g.
#'   `"C"` and `"T"`.
#' @param enzyme_panel `data.frame` with columns `name`, `motif`,
#'   `cut_offset` (default [pira_enzyme_panel()]). You can pass your own
#'   panel with the same columns.
#' @param flank_length How many bases of genomic context to consider on
#'   each side of the SNP when testing motif alignments (default `12`;
#'   should be at least as long as the longest motif in `enzyme_panel`).
#' @param gap Number of unmodifiable genomic bases between a primer's true
#'   3' end and the SNP (default `0`: the SNP is the base immediately
#'   following the primer, the classic PIRA-PCR geometry).
#' @param min_dist_3prime Minimum allowed distance (bp) between the
#'   primer's 3' terminal base and the nearest engineered mismatch, so the
#'   mismatch does not sit at the very last base (which risks blocking
#'   polymerase extension). Default `2`.
#' @param sides Which side(s) of the SNP to search for the engineered
#'   primer: `"upstream"` (a forward primer whose 3' end approaches the
#'   SNP), `"downstream"` (a reverse primer approaching from the other
#'   side), or both (default `c("upstream", "downstream")`).
#'
#' @return A `data.frame`, one row per viable (enzyme, side, frame)
#'   combination, sorted with natural sites first, then by ascending
#'   number of required mismatches, with columns: `enzyme`, `motif`,
#'   `cut_offset` (as given in `enzyme_panel`, i.e. relative to the
#'   enzyme's own canonical orientation), `search_motif`,
#'   `search_cut_offset` (the forward-strand representation to use when
#'   verifying the site on an actual amplicon with
#'   [find_restriction_site()] -- identical to `motif`/`cut_offset` for
#'   `side = "upstream"` or for palindromic motifs, but different for a
#'   non-palindromic enzyme on `side = "downstream"`), `side`,
#'   `snp_motif_pos`, `n_mismatch`, `natural_site`, `cuts_on` (whether the
#'   site forms on the `"ref"` or `"alt"` allele), `mismatch_dist_3prime`,
#'   `required_base`, `original_base` (the last three are comma-separated
#'   strings when `n_mismatch > 1`).
#'
#' @details
#' If any row has `natural_site = TRUE`, a standard (non-PIRA)
#' `design_primers()` + `find_restriction_site()` design should work for
#' that enzyme without engineering anything -- worth checking before
#' committing to a PIRA-PCR primer.
#'
#' This function only identifies *which* (enzyme, mismatch) combinations
#' are chemically viable; it does not yet build full primer pairs. Use
#' [design_pira_primers()] to turn a chosen row into an actual primer pair
#' with Tm/GC/dimer/hairpin filtering and amplicon simulation.
#'
#' @seealso [design_pira_primers()], [pira_enzyme_panel()]
#' @export
find_pira_sites <- function(gene_seq, snp_pos, ref_allele, alt_allele,
                             enzyme_panel = pira_enzyme_panel(),
                             flank_length = 12, gap = 0,
                             min_dist_3prime = 2,
                             sides = c("upstream", "downstream")) {
  seq_char <- if (inherits(gene_seq, "DNAString")) as.character(gene_seq) else gene_seq
  seq_char <- toupper(seq_char)
  ref_allele <- toupper(ref_allele)
  alt_allele <- toupper(alt_allele)

  if (nchar(ref_allele) != 1 || nchar(alt_allele) != 1) {
    stop("ref_allele and alt_allele must each be a single base.")
  }

  n <- nchar(seq_char)
  up_start <- max(1, snp_pos - flank_length)
  down_end <- min(n, snp_pos + flank_length)
  if (snp_pos <= 1 || snp_pos >= n) {
    stop("snp_pos is too close to the edge of gene_seq for the requested flank_length.")
  }

  upstream_fwd   <- strsplit(substr(seq_char, up_start, snp_pos - 1), "")[[1]]
  downstream_fwd <- strsplit(substr(seq_char, snp_pos + 1, down_end), "")[[1]]

  results <- list()

  for (side in sides) {
    if (side == "upstream") {
      up_flank   <- upstream_fwd
      down_flank <- downstream_fwd
      ref_b <- ref_allele
      alt_b <- alt_allele
    } else if (side == "downstream") {
      up_flank   <- .revcomp_vec(downstream_fwd)
      down_flank <- .revcomp_vec(upstream_fwd)
      ref_b <- .complement_base(ref_allele)
      alt_b <- .complement_base(alt_allele)
    } else {
      stop(sprintf("Unrecognized value in 'sides': '%s' (use 'upstream'/'downstream').", side))
    }

    for (i in seq_len(nrow(enzyme_panel))) {
      enz_name <- enzyme_panel$name[i]
      motif_chr <- toupper(enzyme_panel$motif[i])
      motif_vec <- strsplit(motif_chr, "")[[1]]
      L <- length(motif_vec)
      panel_cut_offset <- enzyme_panel$cut_offset[i]

      # For the downstream (mirrored/bottom-strand) side, the physical site
      # reads as the reverse complement of the panel's motif in this frame;
      # using the un-mirrored motif here would silently give wrong results
      # for any non-palindromic enzyme. search_motif/search_cut_offset are
      # the forward-strand representation to use later, when verifying the
      # real cut with find_restriction_site() on the actual amplicon.
      if (side == "upstream") {
        scan_motif_vec <- motif_vec
        search_motif <- motif_chr
        search_cut_offset <- panel_cut_offset
      } else {
        scan_motif_vec <- .iupac_revcomp_motif(motif_vec)
        search_motif <- paste(scan_motif_vec, collapse = "")
        search_cut_offset <- L - panel_cut_offset
      }

      for (s in seq_len(L)) {
        hit <- .pira_scan_frame(
          upstream_flank = up_flank, downstream_flank = down_flank,
          ref_base = ref_b, alt_base = alt_b,
          motif = scan_motif_vec, snp_motif_pos = s,
          min_dist_3prime = min_dist_3prime, gap = gap
        )
        if (is.null(hit)) next

        results[[length(results) + 1]] <- data.frame(
          enzyme = enz_name, motif = motif_chr,
          cut_offset = panel_cut_offset,
          search_motif = search_motif, search_cut_offset = search_cut_offset,
          side = side, snp_motif_pos = s,
          n_mismatch = hit$n_mismatch, natural_site = hit$natural_site,
          cuts_on = hit$cuts_on,
          mismatch_dist_3prime = paste(hit$mismatch_dist_3prime, collapse = ","),
          required_base = paste(hit$required_base, collapse = ","),
          original_base = paste(hit$original_base, collapse = ","),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    message("No viable (enzyme, frame) combination was found with the given panel and constraints. Try increasing flank_length, relaxing min_dist_3prime, or supplying a larger enzyme_panel.")
    return(data.frame(
      enzyme = character(), motif = character(), cut_offset = integer(),
      search_motif = character(), search_cut_offset = integer(),
      side = character(), snp_motif_pos = integer(), n_mismatch = integer(),
      natural_site = logical(), cuts_on = character(),
      mismatch_dist_3prime = character(), required_base = character(),
      original_base = character(), stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, results)
  out <- out[order(-out$natural_site, out$n_mismatch, -nchar(out$mismatch_dist_3prime)), ]
  rownames(out) <- NULL

  n_natural <- sum(out$natural_site)
  if (n_natural > 0) {
    message(sprintf(
      "%d combination(s) form a NATURAL site (0 mismatches) -- ordinary design_primers() + find_restriction_site() may work without PIRA-PCR for those enzymes.",
      n_natural
    ))
  }
  message(sprintf("Found %d viable (enzyme, side, frame) combination(s) in total.", nrow(out)))

  out
}
