# Internal functions used by design_primers(). They are not exported
# because design_primers() is the intended public interface; these
# functions manipulate intermediate structures (raw candidate lists,
# internal column names) that may change without affecting the package's
# public API.

# Generates all possible subsequences (sliding window) of a subregion, for
# lengths between len_min and len_max. Returns a list of lists with the
# subsequence and its local position (1-based, inclusive) within
# dna_sequence.
#' @keywords internal
.extract_candidates <- function(dna_sequence, len_min, len_max) {
  candidates <- vector("list", 0)
  idx <- 1
  for (l in len_min:len_max) {
    windows <- length(dna_sequence) - l + 1
    if (windows > 0) {
      for (i in 1:windows) {
        candidates[[idx]] <- list(
          seq = Biostrings::subseq(dna_sequence, i, i + l - 1),
          local_start = i,
          local_end = i + l - 1
        )
        idx <- idx + 1
      }
    }
  }
  candidates
}

# Filters a list of raw candidates by Tm, %GC, self-dimer and hairpin, and
# returns a data.frame with the absolute genomic position of each
# candidate that passes the filters.
#' @keywords internal
.filter_candidates <- function(candidates_list, region_offset, is_reverse,
                                tm_min, tm_max, gc_min, gc_max,
                                Na, Mg, dNTPs, oligo_conc_nM,
                                dimer_dg_min = -6, hairpin_dg_min = -6) {
  results <- data.frame(
    sequence = character(), base_pairs = integer(), tm = numeric(),
    gc = numeric(), auto_dimer = numeric(), hairpin = numeric(),
    genomic_start = integer(), genomic_end = integer(),
    strand = character(), stringsAsFactors = FALSE
  )

  for (cand in candidates_list) {
    frag <- cand$seq
    g_start <- region_offset + cand$local_start - 1
    g_end <- region_offset + cand$local_end - 1

    if (is_reverse) {
      primer <- Biostrings::reverseComplement(frag)
      strand <- "reverse"
    } else {
      primer <- frag
      strand <- "forward"
    }

    seq_char <- as.character(primer)

    tm <- tryCatch(
      calc_tm(primer, Na = Na, Mg = Mg, dNTPs = dNTPs, oligo_conc_nM = oligo_conc_nM),
      error = function(e) NA_real_
    )
    if (is.na(tm) || tm < tm_min || tm > tm_max) next

    gc <- as.numeric((Biostrings::letterFrequency(primer, "GC") / length(primer)) * 100)
    if (gc < gc_min || gc > gc_max) next

    auto_info <- evaluate_dimer(seq_char, seq_char)
    if (auto_info$dg_approx < dimer_dg_min) next

    hairpin_info <- evaluate_hairpin(seq_char)
    if (hairpin_info$dg_approx < hairpin_dg_min) next

    results <- rbind(results, data.frame(
      sequence = seq_char, base_pairs = length(primer),
      tm = round(tm, 1), gc = round(gc, 1),
      auto_dimer = round(auto_info$dg_approx, 1),
      hairpin = round(hairpin_info$dg_approx, 1),
      genomic_start = g_start, genomic_end = g_end,
      strand = strand, stringsAsFactors = FALSE
    ))
  }

  results
}

# If too many candidates remain per strand, caps the fwd x rev combinatorics
# (which grows as N x M) to a manageable number, keeping the ones with the
# best individual quality (Tm and GC closest to target, lowest
# dimer/hairpin risk).
#' @keywords internal
.preselect_candidates <- function(df, n_max, tm_min, tm_max, gc_min, gc_max) {
  if (nrow(df) <= n_max) return(df)

  tm_target <- mean(c(tm_min, tm_max))
  gc_target <- mean(c(gc_min, gc_max))

  individual_score <- abs(df$tm - tm_target) +
    abs(df$gc - gc_target) * 0.2 +
    pmax(0, -df$auto_dimer - 2) * 0.3 +
    pmax(0, -df$hairpin - 2) * 0.3

  df <- df[order(individual_score), ]
  utils::head(df, n_max)
}

# Generates every forward x reverse combination, filters by amplicon size,
# Tm difference, heterodimer risk, and expected gel resolution after
# digestion (minimum difference between fragments and maximum size of the
# small fragment), and returns a data.frame with the pairs that meet all
# criteria.
#' @keywords internal
.generate_pairs <- function(candidates_f, candidates_r, snp_pos,
                             amplicon_min, amplicon_max, tm_diff_max,
                             heterodimer_dg_min = -6,
                             min_fragment_diff = 25,
                             max_small_fragment = 100,
                             min_fragment_size = 40,
                             verbose = TRUE) {
  pairs <- data.frame()
  fails_amplicon <- 0
  fails_tm <- 0
  fails_dimer <- 0
  fails_rflp <- 0

  for (i in seq_len(nrow(candidates_f))) {
    fwd <- candidates_f[i, ]
    for (j in seq_len(nrow(candidates_r))) {
      rev <- candidates_r[j, ]

      amplicon_len <- rev$genomic_end - fwd$genomic_start + 1

      if (amplicon_len < amplicon_min || amplicon_len > amplicon_max) {
        fails_amplicon <- fails_amplicon + 1
        next
      }

      tm_diff <- abs(fwd$tm - rev$tm)
      if (tm_diff > tm_diff_max) {
        fails_tm <- fails_tm + 1
        next
      }

      hetero <- evaluate_dimer(fwd$sequence, rev$sequence)
      if (hetero$dg_approx < heterodimer_dg_min) {
        fails_dimer <- fails_dimer + 1
        next
      }

      # Expected gel resolution: snp_pos is used as an approximation of the
      # enzyme cut site (valid when the SNP falls within the recognition
      # motif, the typical PCR-RFLP case). The exact cut site can be
      # refined later on the real amplicon with find_restriction_site().
      fragment_1 <- snp_pos - fwd$genomic_start + 1
      fragment_2 <- amplicon_len - fragment_1
      small_fragment <- min(fragment_1, fragment_2)
      fragment_diff <- abs(fragment_1 - fragment_2)

      if (fragment_diff < min_fragment_diff ||
          small_fragment > max_small_fragment ||
          small_fragment < min_fragment_size) {
        fails_rflp <- fails_rflp + 1
        next
      }

      pairs <- rbind(pairs, data.frame(
        forward_seq = fwd$sequence, forward_tm = fwd$tm, forward_gc = fwd$gc,
        forward_start = fwd$genomic_start, forward_end = fwd$genomic_end,
        reverse_seq = rev$sequence, reverse_tm = rev$tm, reverse_gc = rev$gc,
        reverse_start = rev$genomic_start, reverse_end = rev$genomic_end,
        amplicon_bp = amplicon_len, tm_diff = round(tm_diff, 1),
        heterodimer_dg = hetero$dg_approx,
        fragment_1 = fragment_1, fragment_2 = fragment_2,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (verbose) {
    cat("\n--- PAIR FILTERING SUMMARY ---\n")
    cat(sprintf("Failed due to amplicon size:        %d\n", fails_amplicon))
    cat(sprintf("Failed due to Tm difference:         %d\n", fails_tm))
    cat(sprintf("Failed due to (hetero)dimer:         %d\n", fails_dimer))
    cat(sprintf("Failed due to RFLP gel resolution:   %d\n", fails_rflp))
    cat(sprintf("Successful pairs:                    %d\n", nrow(pairs)))
  }

  pairs
}
