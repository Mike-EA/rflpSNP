# Builds the actual engineered-primer sequence (5' -> 3', as it would be
# synthesized) for a chosen (side, frame) PIRA-PCR combination and a given
# primer length, applying the required mismatch(es) at the correct
# position(s) counted from the 3' end. Distances/required bases are taken
# directly from a find_pira_sites() row (they do not depend on the chosen
# primer length, only on the frame).
#' @keywords internal
.build_engineered_primer <- function(gene_seq, snp_pos, side, gap, primer_length,
                                     mismatch_dist_3prime, required_base) {
  seq_char <- if (inherits(gene_seq, "DNAString")) as.character(gene_seq) else gene_seq
  seq_char <- toupper(seq_char)
  n <- nchar(seq_char)

  if (side == "upstream") {
    anchor <- snp_pos - 1 - gap
    span_start <- anchor - primer_length + 1
    span_end <- anchor
    if (span_start < 1) return(NULL)
    primer_bases <- strsplit(substr(seq_char, span_start, span_end), "")[[1]]
  } else if (side == "downstream") {
    anchor <- snp_pos + 1 + gap
    span_start <- anchor
    span_end <- anchor + primer_length - 1
    if (span_end > n) return(NULL)
    fwd_bases <- strsplit(substr(seq_char, span_start, span_end), "")[[1]]
    primer_bases <- .revcomp_vec(fwd_bases)
  } else {
    stop("side must be 'upstream' or 'downstream'")
  }

  L <- length(primer_bases)
  for (k in seq_along(mismatch_dist_3prime)) {
    dist <- mismatch_dist_3prime[k]
    pos <- L - dist + 1
    if (pos < 1 || pos > L) return(NULL)
    primer_bases[pos] <- toupper(required_base[k])
  }

  list(
    sequence = paste(primer_bases, collapse = ""),
    genomic_start = span_start, genomic_end = span_end,
    side = side, length = primer_length
  )
}

# Builds the full amplicon sequence expected for a given allele: the
# region spanned by the engineered primer reflects the primer's actual
# (mismatched) sequence -- since that sequence gets copied into every
# amplicon regardless of genotype once incorporated by the polymerase --
# while the SNP position itself reflects the requested allele, and
# everything else is unmodified genomic sequence.
#' @keywords internal
.build_allele_amplicon <- function(gene_seq, amplicon_start, amplicon_end,
                                   engineered_seq, engineered_span, engineered_side,
                                   snp_pos, allele_base) {
  seq_char <- if (inherits(gene_seq, "DNAString")) as.character(gene_seq) else gene_seq
  seq_char <- toupper(seq_char)
  bases <- strsplit(substr(seq_char, amplicon_start, amplicon_end), "")[[1]]

  eng_start_rel <- engineered_span[1] - amplicon_start + 1
  eng_end_rel   <- engineered_span[2] - amplicon_start + 1

  eng_fwd_seq <- if (engineered_side == "upstream") {
    engineered_seq
  } else {
    paste(.revcomp_vec(strsplit(engineered_seq, "")[[1]]), collapse = "")
  }
  bases[eng_start_rel:eng_end_rel] <- strsplit(eng_fwd_seq, "")[[1]]

  snp_rel <- snp_pos - amplicon_start + 1
  bases[snp_rel] <- toupper(allele_base)

  paste(bases, collapse = "")
}

#' Design a full PIRA-PCR primer pair
#'
#' Builds complete PIRA-PCR (mismatch PCR-RFLP) primer pairs: one primer
#' carries a deliberately introduced mismatch near its 3' end so that,
#' combined with the SNP itself, it creates a restriction site that
#' discriminates the reference and alternate alleles. The other primer is a
#' normal, perfectly-matching primer found with the same candidate
#' generation and filtering machinery used by [design_primers()].
#'
#' @param gene_seq `DNAString` (or character) with the reference sequence.
#' @param snp_pos 1-based coordinate of the SNP in `gene_seq`.
#' @param ref_allele,alt_allele The two alleles at the SNP position
#'   (single-character strings), e.g. `"C"` and `"T"`.
#' @param sites `data.frame` from [find_pira_sites()] describing which
#'   (enzyme, side, frame) combinations to try, in priority order. If
#'   `NULL` (default), it is computed internally by calling
#'   [find_pira_sites()] with `enzyme_panel`, `gap` and `min_dist_3prime`.
#' @param enzyme_panel Enzyme panel used when `sites = NULL` (default
#'   [pira_enzyme_panel()]).
#' @param max_sites_to_try How many of the top rows of `sites` to attempt
#'   building a full primer pair for (default `60`). Not every chemically
#'   viable (enzyme, frame) combination will yield a pair passing the
#'   Tm/GC/amplicon/dimer filters, so trying more than one is recommended.
#' @param gap,min_dist_3prime Passed to [find_pira_sites()] when
#'   `sites = NULL` (see that function).
#' @param search_flank Outer safety cap (bp) on how far to search for the
#'   normal primer on the side opposite the engineered one (default `400`).
#'   The actual search region is narrowed further, and automatically,
#'   based on `amplicon_min`/`amplicon_max` and the engineered primer's own
#'   genomic span -- `search_flank` only guards against searching
#'   arbitrarily far if you set a very large `amplicon_max`. Lower it to
#'   deliberately restrict the search distance for other reasons (e.g. to
#'   avoid a distant repeat region you already know about).
#' @param length_min,length_max Primer length range in bp (default
#'   `18`/`24`).
#' @param tm_min,tm_max,gc_min,gc_max,tm_diff_max Same meaning as in
#'   [design_primers()].
#' @param amplicon_min,amplicon_max Acceptable amplicon size range in bp
#'   (default `150`/`300`).
#' @param min_fragment_diff,max_small_fragment,min_fragment_size Same
#'   meaning as in [design_primers()], but applied to the fragment sizes
#'   actually observed after simulating digestion of both allele amplicons
#'   (see Details). In PIRA-PCR the cut site sits immediately next to the
#'   engineered primer's 3' end, so one fragment is systematically small
#'   (roughly the engineered primer's length plus the enzyme's cut offset)
#'   -- this is expected for the technique, not a design flaw, and the
#'   achievable size is often dictated by which primer lengths satisfy
#'   `tm_min`/`tm_max` rather than being freely choosable. For that reason
#'   `min_fragment_size` here is a true detectability floor (default `10`;
#'   below this a fragment is not realistically callable on any gel), not a
#'   comfort threshold -- designs with a small fragment above the floor but
#'   still on the small side are not rejected, only ranked lower via `score`
#'   (see `small_fragment` in the returned `top`/`best_pair`). Raise
#'   `min_fragment_size` if your gel setup needs a larger minimum.
#' @param Na,Mg,dNTPs,oligo_conc_nM PCR reaction parameters used for Tm
#'   (see [calc_tm()]). Note the Tm of the engineered primer is
#'   necessarily approximate, since it is computed as if the primer were a
#'   perfect match to itself; the true duplex against the genomic template
#'   contains an intentional mismatch. Verify with a dedicated tool (e.g.
#'   IDT OligoAnalyzer, entering the mismatch explicitly) before ordering
#'   synthesis.
#' @param n_top Number of best pairs to report (default `5`).
#' @param verbose If `TRUE` (default), prints progress messages.
#'
#' @return An object of class `"rflp_pira_primers"` with `top`, `best_pair`,
#'   `n_pairs_total` and `sites_tried`, analogous to the object returned by
#'   [design_primers()], plus PIRA-specific columns: `enzyme`, `side`,
#'   `n_mismatch`, `cuts_on`, `n_incidental_sites` (extra, SNP-unrelated
#'   cut sites elsewhere in the amplicon -- present regardless of genotype,
#'   so they don't confound genotyping but add extra bands), `small_fragment`
#'   (the smaller of the two digestion fragments), `engineered_primer`,
#'   `mismatch_positions_from_3prime`, `required_base`, `original_base`.
#'
#' @details
#' For every candidate pair, both allele amplicons are built explicitly
#' (substituting the engineered primer's own sequence into its genomic
#' span, and the requested allele at the SNP position) and passed through
#' [find_restriction_site()] for real -- the cut site and fragment sizes
#' are verified, not just assumed from the frame geometry. A pair is
#' rejected if the site does not appear on the intended allele, or if it
#' unexpectedly appears on both (which would mean it no longer
#' discriminates the two alleles).
#'
#' @seealso [find_pira_sites()], [design_primers()], [simulate_pcr()],
#'   [find_restriction_site()]
#' @export
design_pira_primers <- function(gene_seq, snp_pos, ref_allele, alt_allele,
                                sites = NULL, enzyme_panel = pira_enzyme_panel(),
                                max_sites_to_try = 60,
                                gap = 0, min_dist_3prime = 2,
                                search_flank = 400,
                                length_min = 18, length_max = 24,
                                tm_min = 50, tm_max = 65,
                                gc_min = 35, gc_max = 65,
                                tm_diff_max = 5,
                                amplicon_min = 150, amplicon_max = 300,
                                min_fragment_diff = 25, max_small_fragment = 100,
                                min_fragment_size = 10,
                                Na = 100, Mg = 2, dNTPs = 0.2, oligo_conc_nM = 500,
                                n_top = 5, verbose = TRUE) {

  ref_allele <- toupper(ref_allele)
  alt_allele <- toupper(alt_allele)

  if (is.null(sites)) {
    sites <- find_pira_sites(
      gene_seq, snp_pos, ref_allele, alt_allele,
      enzyme_panel = enzyme_panel,
      flank_length = max(length_max + gap + 2, 12),
      gap = gap, min_dist_3prime = min_dist_3prime
    )
  }
  if (nrow(sites) == 0) {
    stop("No viable PIRA-PCR (enzyme, frame) combination was found. Try a larger/custom enzyme_panel, a smaller min_dist_3prime, or inspect find_pira_sites() directly.")
  }
  sites_to_try <- utils::head(sites, max_sites_to_try)

  seq_char <- if (inherits(gene_seq, "DNAString")) as.character(gene_seq) else gene_seq
  seq_char <- toupper(seq_char)
  n_total <- nchar(seq_char)

  gc_content <- function(s) 100 * nchar(gsub("[^GC]", "", s)) / nchar(s)

  all_pairs <- list()

  # Diagnostic counters, aggregated across every site attempted -- printed
  # at the end (always, not just when verbose) if no pair is found, so a
  # failure is never a silent dead end.
  diag <- list(
    sites_tried = 0L,
    sites_no_engineered_primer = 0L,
    sites_no_normal_region = 0L,
    sites_no_normal_primer = 0L,
    pairs_considered = 0L,
    fail_amplicon_size = 0L,
    fail_tm_diff = 0L,
    fail_heterodimer = 0L,
    fail_site_not_on_target_allele = 0L,
    fail_site_on_both_alleles = 0L,
    fail_cut_outside_amplicon = 0L,
    fail_fragment_resolution = 0L
  )

  for (i in seq_len(nrow(sites_to_try))) {
    site <- sites_to_try[i, ]
    diag$sites_tried <- diag$sites_tried + 1L
    mismatch_dist <- if (nchar(site$mismatch_dist_3prime) == 0) {
      integer(0)
    } else as.integer(strsplit(site$mismatch_dist_3prime, ",")[[1]])
    required_base <- if (nchar(site$required_base) == 0) {
      character(0)
    } else strsplit(site$required_base, ",")[[1]]

    if (verbose) cat(sprintf("Trying site %d/%d: %s (%s side, %d mismatch)...\n", i, nrow(sites_to_try), site$enzyme, site$side, site$n_mismatch))

    # 1) engineered-primer candidates over the requested length range
    eng_rows <- list()
    for (Len in length_min:length_max) {
      built <- .build_engineered_primer(seq_char, snp_pos, site$side, gap, Len, mismatch_dist, required_base)
      if (is.null(built)) next
      tm <- tryCatch(calc_tm(built$sequence, Na = Na, Mg = Mg, dNTPs = dNTPs, oligo_conc_nM = oligo_conc_nM), error = function(e) NA_real_)
      if (is.na(tm) || tm < tm_min || tm > tm_max) next
      gc <- gc_content(built$sequence)
      if (gc < gc_min || gc > gc_max) next
      auto_info <- evaluate_dimer(built$sequence, built$sequence)
      if (auto_info$dg_approx < -6) next
      hairpin_info <- evaluate_hairpin(built$sequence)
      if (hairpin_info$dg_approx < -6) next

      eng_rows[[length(eng_rows) + 1]] <- data.frame(
        sequence = built$sequence, base_pairs = Len, tm = round(tm, 1), gc = round(gc, 1),
        genomic_start = built$genomic_start, genomic_end = built$genomic_end,
        stringsAsFactors = FALSE
      )
    }
    if (length(eng_rows) == 0) {
      diag$sites_no_engineered_primer <- diag$sites_no_engineered_primer + 1L
      next
    }
    eng_df <- do.call(rbind, eng_rows)
    eng_df <- utils::head(eng_df[order(abs(eng_df$tm - mean(c(tm_min, tm_max)))), ], 15)

    # 2) normal-primer candidates on the opposite side, reusing existing
    # infra. Narrow the search region to genomic positions that could
    # actually yield an amplicon within [amplicon_min, amplicon_max],
    # given the engineered primer's real (already Tm/GC-filtered) anchor
    # range in eng_df. Searching blindly over a flat search_flank window
    # wastes most of the work on positions that can never produce a valid
    # amplicon size (too close to the SNP) while sometimes missing
    # positions needed for longer engineered primers (too far) if
    # search_flank isn't generous enough.
    if (site$side == "upstream") {
      eng_start_min <- min(eng_df$genomic_start)
      eng_start_max <- max(eng_df$genomic_start)
      valid_end_lo <- eng_start_min + amplicon_min - 1
      valid_end_hi <- eng_start_max + amplicon_max - 1

      norm_start <- max(snp_pos + 1 + gap, valid_end_lo - length_max + 1)
      norm_end   <- min(n_total, snp_pos + search_flank, valid_end_hi)
    } else {
      eng_end_min <- min(eng_df$genomic_end)
      eng_end_max <- max(eng_df$genomic_end)
      valid_start_lo <- eng_end_min - amplicon_max + 1
      valid_start_hi <- eng_end_max - amplicon_min + 1

      norm_end   <- min(snp_pos - 1 - gap, valid_start_hi + length_max - 1)
      norm_start <- max(1, snp_pos - search_flank, valid_start_lo)
    }
    if (norm_end - norm_start + 1 < length_min) {
      diag$sites_no_normal_region <- diag$sites_no_normal_region + 1L
      next
    }

    norm_region <- Biostrings::DNAString(substr(seq_char, norm_start, norm_end))
    norm_raw <- .extract_candidates(norm_region, length_min, length_max)
    norm_df <- .filter_candidates(
      norm_raw, region_offset = norm_start, is_reverse = (site$side == "upstream"),
      tm_min = tm_min, tm_max = tm_max, gc_min = gc_min, gc_max = gc_max,
      Na = Na, Mg = Mg, dNTPs = dNTPs, oligo_conc_nM = oligo_conc_nM
    )
    if (nrow(norm_df) == 0) {
      diag$sites_no_normal_primer <- diag$sites_no_normal_primer + 1L
      next
    }
    # Preselect toward the engineered primers' actual Tm, not the generic
    # midpoint of [tm_min, tm_max] -- otherwise we systematically keep
    # normal candidates whose Tm is close to the *range's* center even when
    # every viable engineered primer's Tm sits well off to one side,
    # guaranteeing most pairs fail the tm_diff_max check downstream.
    eng_tm_center <- mean(eng_df$tm)
    norm_df <- .preselect_candidates(norm_df, 40, eng_tm_center, eng_tm_center, gc_min, gc_max)

    # 3) pair engineered x normal, verifying the real digestion outcome
    for (a in seq_len(nrow(eng_df))) {
      eng <- eng_df[a, ]
      for (b in seq_len(nrow(norm_df))) {
        norm <- norm_df[b, ]
        diag$pairs_considered <- diag$pairs_considered + 1L

        if (site$side == "upstream") {
          amp_start <- eng$genomic_start; amp_end <- norm$genomic_end
          fwd_seq <- eng$sequence; fwd_tm <- eng$tm; fwd_gc <- eng$gc
          rev_seq <- norm$sequence; rev_tm <- norm$tm; rev_gc <- norm$gc
        } else {
          amp_start <- norm$genomic_start; amp_end <- eng$genomic_end
          fwd_seq <- norm$sequence; fwd_tm <- norm$tm; fwd_gc <- norm$gc
          rev_seq <- eng$sequence; rev_tm <- eng$tm; rev_gc <- eng$gc
        }
        amplicon_len <- amp_end - amp_start + 1
        if (amplicon_len < amplicon_min || amplicon_len > amplicon_max) {
          diag$fail_amplicon_size <- diag$fail_amplicon_size + 1L
          next
        }

        tm_diff <- abs(eng$tm - norm$tm)
        if (tm_diff > tm_diff_max) {
          diag$fail_tm_diff <- diag$fail_tm_diff + 1L
          next
        }

        hetero <- evaluate_dimer(eng$sequence, norm$sequence)
        if (hetero$dg_approx < -6) {
          diag$fail_heterodimer <- diag$fail_heterodimer + 1L
          next
        }

        amp_ref_seq <- .build_allele_amplicon(
          seq_char, amp_start, amp_end,
          engineered_seq = eng$sequence,
          engineered_span = c(eng$genomic_start, eng$genomic_end),
          engineered_side = site$side, snp_pos = snp_pos, allele_base = ref_allele
        )
        amp_alt_seq <- .build_allele_amplicon(
          seq_char, amp_start, amp_end,
          engineered_seq = eng$sequence,
          engineered_span = c(eng$genomic_start, eng$genomic_end),
          engineered_side = site$side, snp_pos = snp_pos, allele_base = alt_allele
        )

        # Use search_motif/search_cut_offset (the forward-strand
        # representation), not the panel's raw motif/cut_offset -- for a
        # "downstream" site these differ whenever the enzyme is not
        # palindromic.
        site_ref <- find_restriction_site(Biostrings::DNAString(amp_ref_seq), enzyme_motif = site$search_motif, cut_offset = site$search_cut_offset)
        site_alt <- find_restriction_site(Biostrings::DNAString(amp_alt_seq), enzyme_motif = site$search_motif, cut_offset = site$search_cut_offset)

        cutting <- if (site$cuts_on == "ref") site_ref else site_alt
        noncutting <- if (site$cuts_on == "ref") site_alt else site_ref

        # find_restriction_site() returns positions local to the amplicon
        # string (1..amplicon_len), not genomic coordinates -- amp_start
        # must NOT be subtracted from them.
        snp_rel <- snp_pos - amp_start + 1

        if (cutting$n_sites == 0) {
          diag$fail_site_not_on_target_allele <- diag$fail_site_not_on_target_allele + 1L
          next
        }

        # Identify the engineered/diagnostic occurrence specifically: the
        # one whose span covers the SNP position. Short motifs (e.g. 4bp)
        # commonly occur elsewhere in the amplicon by chance; such
        # incidental sites are present on BOTH alleles (they are unrelated
        # to the SNP) and must not, by themselves, disqualify the pair.
        cutting_at_snp <- cutting$sites[cutting$sites$start <= snp_rel & cutting$sites$end >= snp_rel, , drop = FALSE]
        if (nrow(cutting_at_snp) == 0) {
          diag$fail_site_not_on_target_allele <- diag$fail_site_not_on_target_allele + 1L
          next
        }

        noncutting_at_snp <- noncutting$sites[noncutting$sites$start <= snp_rel & noncutting$sites$end >= snp_rel, , drop = FALSE]
        if (nrow(noncutting_at_snp) > 0) {
          diag$fail_site_on_both_alleles <- diag$fail_site_on_both_alleles + 1L
          next
        }

        n_incidental_sites <- (cutting$n_sites - 1) + noncutting$n_sites

        cut_pos <- cutting_at_snp$cut_pos[1]
        if (cut_pos < 1 || cut_pos >= amplicon_len) {
          diag$fail_cut_outside_amplicon <- diag$fail_cut_outside_amplicon + 1L
          next
        }

        fragment_1 <- cut_pos
        fragment_2 <- amplicon_len - fragment_1
        small_fragment <- min(fragment_1, fragment_2)
        if (abs(fragment_1 - fragment_2) < min_fragment_diff ||
            small_fragment > max_small_fragment ||
            small_fragment < min_fragment_size) {
          diag$fail_fragment_resolution <- diag$fail_fragment_resolution + 1L
          next
        }

        all_pairs[[length(all_pairs) + 1]] <- data.frame(
          enzyme = site$enzyme, side = site$side, n_mismatch = site$n_mismatch,
          cuts_on = site$cuts_on, n_incidental_sites = n_incidental_sites,
          forward_seq = fwd_seq, forward_tm = fwd_tm, forward_gc = fwd_gc,
          reverse_seq = rev_seq, reverse_tm = rev_tm, reverse_gc = rev_gc,
          amplicon_bp = amplicon_len, tm_diff = round(tm_diff, 1),
          heterodimer_dg = hetero$dg_approx,
          fragment_1 = fragment_1, fragment_2 = fragment_2, small_fragment = small_fragment,
          engineered_primer = eng$sequence,
          mismatch_positions_from_3prime = site$mismatch_dist_3prime,
          required_base = site$required_base, original_base = site$original_base,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(all_pairs) == 0) {
    cat("\n--- PIRA-PCR DESIGN FAILURE DIAGNOSTIC ---\n")
    cat(sprintf("Sites attempted:                                  %d (of %d available)\n", diag$sites_tried, nrow(sites)))
    cat(sprintf("  ...with no engineered primer in Tm/GC/dimer range: %d\n", diag$sites_no_engineered_primer))
    cat(sprintf("  ...with no room for a normal primer on the other side: %d\n", diag$sites_no_normal_region))
    cat(sprintf("  ...with no normal primer passing Tm/GC filtering: %d\n", diag$sites_no_normal_primer))
    cat(sprintf("Engineered x normal pairs considered:             %d\n", diag$pairs_considered))
    cat(sprintf("  ...rejected: amplicon size out of range:        %d\n", diag$fail_amplicon_size))
    cat(sprintf("  ...rejected: Tm difference too large:           %d\n", diag$fail_tm_diff))
    cat(sprintf("  ...rejected: heterodimer risk:                  %d\n", diag$fail_heterodimer))
    cat(sprintf("  ...rejected: site absent on target allele:      %d\n", diag$fail_site_not_on_target_allele))
    cat(sprintf("  ...rejected: site present on BOTH alleles:      %d\n", diag$fail_site_on_both_alleles))
    cat(sprintf("  ...rejected: cut position outside amplicon:     %d\n", diag$fail_cut_outside_amplicon))
    cat(sprintf("  ...rejected: fragment sizes not resolvable:     %d\n", diag$fail_fragment_resolution))
    cat("-------------------------------------------\n")

    if (diag$sites_no_engineered_primer == diag$sites_tried) {
      cat("Every attempted site failed at the engineered-primer stage: no length between length_min and length_max produced a primer within [tm_min, tm_max] / [gc_min, gc_max] / dimer / hairpin thresholds. This is the classic PIRA-PCR bottleneck (the primer's 3' end is anchored, so only its length can vary). Try: a larger enzyme_panel (more independent frames to attempt), widening tm_min/tm_max or gc_min/gc_max, or increasing max_sites_to_try.\n")
    } else if (diag$fail_fragment_resolution > 0 && nrow(sites) > 0 && diag$pairs_considered == diag$fail_fragment_resolution) {
      cat("Pairs were found and the engineered site verified correctly, but the two fragments were not distinguishable enough. Try relaxing min_fragment_diff / max_small_fragment / min_fragment_size.\n")
    }

    stop("No valid PIRA-PCR primer pair was found for the top candidate sites. See the diagnostic above for where it failed.")
  }

  result_df <- do.call(rbind, all_pairs)
  tm_target <- mean(c(tm_min, tm_max))
  gc_target <- mean(c(gc_min, gc_max))
  result_df$score <- with(result_df,
                          abs(forward_tm - tm_target) + abs(reverse_tm - tm_target) +
                            abs(forward_gc - gc_target) * 0.2 + abs(reverse_gc - gc_target) * 0.2 +
                            tm_diff * 1.5 + pmax(0, -heterodimer_dg - 4) * 1.0 + n_mismatch * 2.0 +
                            n_incidental_sites * 0.5 + pmax(0, 25 - small_fragment) * 0.3
  )
  result_df <- result_df[order(result_df$score), ]
  top <- utils::head(result_df, n_top)
  top$recommended <- c("YES", rep("NO", nrow(top) - 1))
  rownames(top) <- NULL

  out <- list(
    top = top, best_pair = top[1, ],
    n_pairs_total = nrow(result_df), sites_tried = sites_to_try
  )
  class(out) <- "rflp_pira_primers"

  if (verbose) {
    cat("\n=== Best PIRA-PCR primer pair selected ===\n")
    bp <- out$best_pair
    cat(sprintf("Enzyme: %s | side engineered: %s | mismatches: %d | cuts on: %s allele | incidental sites: %d\n",
                bp$enzyme, bp$side, bp$n_mismatch, bp$cuts_on, bp$n_incidental_sites))
    cat(sprintf("Forward: %s (Tm=%.1f\u00B0C, GC=%.1f%%)\n", bp$forward_seq, bp$forward_tm, bp$forward_gc))
    cat(sprintf("Reverse: %s (Tm=%.1f\u00B0C, GC=%.1f%%)\n", bp$reverse_seq, bp$reverse_tm, bp$reverse_gc))
    cat(sprintf("Expected amplicon: %d bp | Fragments: %d + %d bp\n", bp$amplicon_bp, bp$fragment_1, bp$fragment_2))
  }

  out
}

#' @export
print.rflp_pira_primers <- function(x, ...) {
  cat("<rflp_pira_primers>\n")
  cat(sprintf("%d total candidate pair(s) evaluated; showing the best %d.\n", x$n_pairs_total, nrow(x$top)))
  print(x$top[, c(
    "enzyme", "side", "n_mismatch", "n_incidental_sites", "forward_seq", "reverse_seq",
    "amplicon_bp", "fragment_1", "fragment_2", "small_fragment", "score", "recommended"
  )], row.names = FALSE)
  invisible(x)
}
