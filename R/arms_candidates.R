# Candidate generation, filtering and ranking for tetra-primer ARMS-PCR.

#' @keywords internal
.arms_validate_inputs <- function(gene_seq, snp_pos, ref_allele, alt_allele) {
  if (!inherits(gene_seq, "DNAString")) {
    stop("'gene_seq' must be a Biostrings::DNAString.")
  }
  if (length(snp_pos) != 1L || is.na(snp_pos) || snp_pos != as.integer(snp_pos) ||
      snp_pos < 1L || snp_pos > length(gene_seq)) {
    stop("'snp_pos' must be one 1-based coordinate within 'gene_seq'.")
  }
  ref_allele <- toupper(ref_allele)
  alt_allele <- toupper(alt_allele)
  if (!grepl("^[ACGT]$", ref_allele) || !grepl("^[ACGT]$", alt_allele) ||
      identical(ref_allele, alt_allele)) {
    stop("'ref_allele' and 'alt_allele' must be different single bases (A, C, G or T).")
  }
  observed <- as.character(gene_seq[snp_pos])
  if (!identical(observed, ref_allele)) {
    stop("'ref_allele' does not match the base at 'snp_pos' in 'gene_seq'.")
  }
  list(snp_pos = as.integer(snp_pos), ref_allele = ref_allele, alt_allele = alt_allele)
}

#' @keywords internal
.arms_replace_base <- function(sequence, index, base) {
  paste0(substr(sequence, 1L, index - 1L), base, substr(sequence, index + 1L, nchar(sequence)))
}

#' @keywords internal
.arms_empty_candidates <- function() {
  data.frame(
    sequence = character(), genomic_start = integer(), genomic_end = integer(),
    strand = character(), primer_kind = character(), allele = character(),
    mismatch_3prime = character(), mismatch_index = integer(),
    mismatch_from = character(), mismatch_to = character(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.arms_outer_candidates <- function(gene_seq, starts, ends, strand) {
  rows <- lapply(seq_along(starts), function(i) {
    template <- Biostrings::subseq(gene_seq, starts[i], ends[i])
    sequence <- if (identical(strand, "forward")) template else Biostrings::reverseComplement(template)
    data.frame(
      sequence = as.character(sequence), genomic_start = starts[i], genomic_end = ends[i],
      strand = strand, primer_kind = "outer", allele = NA_character_,
      mismatch_3prime = NA_character_, mismatch_index = NA_integer_,
      mismatch_from = NA_character_, mismatch_to = NA_character_, stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) return(.arms_empty_candidates())
  do.call(rbind, rows)
}

#' @keywords internal
.arms_inner_candidates <- function(gene_seq, snp_pos, lengths, strand, allele,
                                   mismatch_positions) {
  rows <- vector("list", 0L)
  index <- 1L
  terminal <- as.character(Biostrings::complement(Biostrings::DNAString(allele)))
  if (identical(strand, "forward")) terminal <- allele

  for (primer_length in lengths) {
    start <- if (identical(strand, "forward")) snp_pos - primer_length + 1L else snp_pos
    end <- if (identical(strand, "forward")) snp_pos else snp_pos + primer_length - 1L
    if (start < 1L || end > length(gene_seq)) next

    template <- Biostrings::subseq(gene_seq, start, end)
    sequence <- as.character(if (identical(strand, "forward")) template else Biostrings::reverseComplement(template))
    sequence <- .arms_replace_base(sequence, nchar(sequence), terminal)

    for (distance in mismatch_positions) {
      mismatch_index <- nchar(sequence) - distance
      if (mismatch_index < 1L) next
      original <- substr(sequence, mismatch_index, mismatch_index)
      for (replacement in setdiff(c("A", "C", "G", "T"), original)) {
        rows[[index]] <- data.frame(
          sequence = .arms_replace_base(sequence, mismatch_index, replacement),
          genomic_start = start, genomic_end = end, strand = strand,
          primer_kind = "inner", allele = allele,
          mismatch_3prime = paste0("3'-", distance), mismatch_index = mismatch_index,
          mismatch_from = original, mismatch_to = replacement,
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
  }
  if (!length(rows)) return(.arms_empty_candidates())
  do.call(rbind, rows)
}

#' Generate unfiltered tetra-primer ARMS-PCR candidate pools
#'
#' Creates geometrically valid outer and allele-specific inner primer pools.
#' It does not evaluate Tm, GC, secondary structures, dimers, amplicon sizes
#' or ranking; those decisions belong to the stage-2 design function.
#'
#' @noRd
#' @keywords internal
.generate_arms_candidate_pools <- function(gene_seq, snp_pos, ref_allele, alt_allele,
                                            outer_flank = 160L,
                                            outer_length_min = 18L,
                                            outer_length_max = 24L,
                                            inner_length_min = 18L,
                                            inner_length_max = 24L,
                                            mismatch_positions = c(2L, 3L)) {
  input <- .arms_validate_inputs(gene_seq, snp_pos, ref_allele, alt_allele)
  if (outer_flank < 1L || outer_length_min < 1L || outer_length_max < outer_length_min ||
      inner_length_min < 1L || inner_length_max < inner_length_min ||
      !all(mismatch_positions %in% c(2L, 3L))) {
    stop("Invalid candidate-length, flank or mismatch-position parameters.")
  }
  snp_pos <- input$snp_pos
  upstream_start <- max(1L, snp_pos - as.integer(outer_flank))
  downstream_end <- min(length(gene_seq), snp_pos + as.integer(outer_flank))

  outer_fwd_starts <- outer_fwd_ends <- integer()
  outer_rev_starts <- outer_rev_ends <- integer()
  for (primer_length in seq.int(outer_length_min, outer_length_max)) {
    fwd_first_end <- upstream_start + primer_length - 1L
    fwd_end <- if (fwd_first_end <= snp_pos - 1L) seq.int(fwd_first_end, snp_pos - 1L) else integer()
    rev_last_start <- downstream_end - primer_length + 1L
    rev_start <- if (snp_pos + 1L <= rev_last_start) seq.int(snp_pos + 1L, rev_last_start) else integer()
    if (length(fwd_end)) {
      outer_fwd_starts <- c(outer_fwd_starts, fwd_end - primer_length + 1L)
      outer_fwd_ends <- c(outer_fwd_ends, fwd_end)
    }
    if (length(rev_start)) {
      outer_rev_starts <- c(outer_rev_starts, rev_start)
      outer_rev_ends <- c(outer_rev_ends, rev_start + primer_length - 1L)
    }
  }

  pools <- list(
    outer_forward = .arms_outer_candidates(gene_seq, outer_fwd_starts, outer_fwd_ends, "forward"),
    outer_reverse = .arms_outer_candidates(gene_seq, outer_rev_starts, outer_rev_ends, "reverse"),
    ref_forward = .arms_inner_candidates(gene_seq, snp_pos, seq.int(inner_length_min, inner_length_max), "forward", input$ref_allele, mismatch_positions),
    ref_reverse = .arms_inner_candidates(gene_seq, snp_pos, seq.int(inner_length_min, inner_length_max), "reverse", input$ref_allele, mismatch_positions),
    alt_forward = .arms_inner_candidates(gene_seq, snp_pos, seq.int(inner_length_min, inner_length_max), "forward", input$alt_allele, mismatch_positions),
    alt_reverse = .arms_inner_candidates(gene_seq, snp_pos, seq.int(inner_length_min, inner_length_max), "reverse", input$alt_allele, mismatch_positions)
  )
  layouts <- data.frame(
    layout = c("ref_forward_alt_reverse", "alt_forward_ref_reverse"),
    ref_pool = c("ref_forward", "ref_reverse"),
    alt_pool = c("alt_reverse", "alt_forward"),
    stringsAsFactors = FALSE
  )
  set_counts <- vapply(seq_len(nrow(layouts)), function(i) {
    prod(as.numeric(c(nrow(pools$outer_forward), nrow(pools$outer_reverse),
                      nrow(pools[[layouts$ref_pool[i]]]), nrow(pools[[layouts$alt_pool[i]]]))))
  }, numeric(1))
  diagnostics <- data.frame(
    pool = names(pools),
    n_candidates = vapply(pools, nrow, integer(1)),
    stringsAsFactors = FALSE
  )
  result <- list(
    top = data.frame(),
    best_set = NULL,
    candidate_pools = pools,
    layouts = layouts,
    n_geometric_sets = sum(set_counts),
    n_geometric_sets_by_layout = stats::setNames(set_counts, layouts$layout),
    parameters = c(input, list(
      outer_flank = outer_flank,
      outer_length_min = outer_length_min,
      outer_length_max = outer_length_max,
      inner_length_min = inner_length_min,
      inner_length_max = inner_length_max,
      mismatch_positions = mismatch_positions
    )),
    diagnostics = diagnostics
  )
  class(result) <- "rflp_arms_primers"
  result
}

#' @export
print.rflp_arms_primers <- function(x, ...) {
  cat("<rflp_arms_primers>\n")
  if (!is.null(x$n_valid_sets)) {
    cat(sprintf("Geometric sets: %s | valid sets: %s\n",
                format(x$n_geometric_sets, big.mark = ","),
                format(x$n_valid_sets, big.mark = ",")))
    if (x$n_valid_sets) {
      print(x$top[, c("layout", "ref_amplicon_bp", "alt_amplicon_bp",
                     "control_amplicon_bp", "score", "recommended")], row.names = FALSE)
    } else {
      cat("No set passed stage-2 filtering; inspect $diagnostics for exclusions.\n")
    }
  } else if (!x$n_geometric_sets) {
    cat("No complete four-primer set can be formed with the supplied geometry.\n")
  } else {
    cat("Filtering and ranking are pending stage 2.\n")
  }
  invisible(x)
}

#' @keywords internal
.arms_gc <- function(sequence) {
  100 * sum(strsplit(sequence, "", fixed = TRUE)[[1]] %in% c("G", "C")) / nchar(sequence)
}

#' @keywords internal
.arms_evaluate_pool <- function(pool, tm_min, tm_max, gc_min, gc_max,
                                 dimer_dg_min, hairpin_dg_min,
                                 Na, Mg, dNTPs, oligo_conc_nM,
                                 max_raw_candidates = Inf) {
  reasons <- c(tm = 0L, gc = 0L, self_dimer = 0L, hairpin = 0L)
  if (nrow(pool) > max_raw_candidates) {
    pool <- pool[unique(round(seq(1, nrow(pool), length.out = max_raw_candidates))), , drop = FALSE]
  }
  if (!nrow(pool)) return(list(candidates = pool, reasons = reasons))
  kept <- vector("list", 0L)
  for (i in seq_len(nrow(pool))) {
    candidate <- pool[i, , drop = FALSE]
    tm <- tryCatch(calc_tm(candidate$sequence, Na = Na, Mg = Mg, dNTPs = dNTPs,
                           oligo_conc_nM = oligo_conc_nM), error = function(e) NA_real_)
    if (is.na(tm) || tm < tm_min || tm > tm_max) { reasons[["tm"]] <- reasons[["tm"]] + 1L; next }
    gc <- .arms_gc(candidate$sequence)
    if (gc < gc_min || gc > gc_max) { reasons[["gc"]] <- reasons[["gc"]] + 1L; next }
    self_dimer <- if (is.infinite(dimer_dg_min) && dimer_dg_min < 0) 0 else
      evaluate_dimer(candidate$sequence, candidate$sequence)$dg_approx
    if (self_dimer < dimer_dg_min) { reasons[["self_dimer"]] <- reasons[["self_dimer"]] + 1L; next }
    hairpin <- if (is.infinite(hairpin_dg_min) && hairpin_dg_min < 0) 0 else
      evaluate_hairpin(candidate$sequence)$dg_approx
    if (hairpin < hairpin_dg_min) { reasons[["hairpin"]] <- reasons[["hairpin"]] + 1L; next }
    candidate$length_bp <- nchar(candidate$sequence)
    candidate$tm <- round(tm, 1)
    candidate$gc <- round(gc, 1)
    candidate$self_dimer_dg <- self_dimer
    candidate$hairpin_dg <- hairpin
    kept[[length(kept) + 1L]] <- candidate
  }
  candidates <- if (length(kept)) do.call(rbind, kept) else pool[FALSE, ]
  list(candidates = candidates, reasons = reasons)
}

#' @keywords internal
.arms_preselect_pool <- function(pool, n_max, tm_target, gc_target) {
  if (nrow(pool) <= n_max) return(pool)
  score <- abs(pool$tm - tm_target) + abs(pool$gc - gc_target) * 0.2 +
    pmax(0, -pool$self_dimer_dg - 2) * 0.3 + pmax(0, -pool$hairpin_dg - 2) * 0.3
  pool[order(score, pool$sequence), , drop = FALSE][seq_len(n_max), , drop = FALSE]
}

#' @keywords internal
.arms_inner_is_specific <- function(primer, gene_seq, snp_pos, target_allele, nontarget_allele) {
  template <- as.character(Biostrings::subseq(gene_seq, primer$genomic_start, primer$genomic_end))
  snp_index <- snp_pos - primer$genomic_start + 1L
  target_template <- .arms_replace_base(template, snp_index, target_allele)
  nontarget_template <- .arms_replace_base(template, snp_index, nontarget_allele)
  binding_sequence <- function(x) if (identical(primer$strand, "forward")) x else
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  target_binding <- binding_sequence(target_template)
  nontarget_binding <- binding_sequence(nontarget_template)
  if (nchar(primer$sequence) != nchar(target_binding) ||
      nchar(primer$sequence) != nchar(nontarget_binding)) return(FALSE)
  terminal_ok <- identical(substr(primer$sequence, nchar(primer$sequence), nchar(primer$sequence)),
                           substr(target_binding, nchar(target_binding), nchar(target_binding)))
  deliberate_ok <- substr(primer$sequence, primer$mismatch_index, primer$mismatch_index) !=
    substr(target_binding, primer$mismatch_index, primer$mismatch_index)
  target_mismatches <- sum(strsplit(primer$sequence, "", fixed = TRUE)[[1]] !=
                             strsplit(target_binding, "", fixed = TRUE)[[1]])
  nontarget_mismatches <- sum(strsplit(primer$sequence, "", fixed = TRUE)[[1]] !=
                                strsplit(nontarget_binding, "", fixed = TRUE)[[1]])
  terminal_ok && deliberate_ok && target_mismatches == 1L && nontarget_mismatches == 2L
}

#' Design and rank tetra-primer ARMS-PCR assays
#'
#' Builds candidate pools, filters each primer by physicochemical properties,
#' verifies the intended two-mismatch 3' discrimination, and ranks complete
#' four-primer sets.  A result with zero valid sets is returned with exclusion
#' diagnostics so unsuitable SNP contexts remain explainable.
#'
#' @param gene_seq Reference `Biostrings::DNAString` in 5' to 3' orientation.
#' @param snp_pos One-based SNP coordinate.
#' @param ref_allele,alt_allele Reference and explicit alternative base.
#' @param outer_flank Search distance for outer primers.
#' @param outer_length_min,outer_length_max,inner_length_min,inner_length_max Primer length limits.
#' @param mismatch_positions Deliberate mismatch positions from the 3' end.
#' @param tm_min,tm_max,gc_min,gc_max Individual primer limits.
#' @param dimer_dg_min,hairpin_dg_min Heuristic individual structure limits.
#' @param heterodimer_dg_min Heuristic cross-dimer limit for every primer pair.
#' @param control_amplicon_min,control_amplicon_max,allele_amplicon_min,allele_amplicon_max Product size limits.
#' @param min_band_diff Minimum pairwise separation among the three products.
#' @param max_candidates_per_pool Cap after individual filtering to bound combinatorics.
#' @param max_raw_candidates_per_pool Optional cap before physicochemical
#'   evaluation, useful for fast exploratory searches in very large windows.
#' @param n_top Number of highest-ranked sets returned.
#' @param Na,Mg,dNTPs,oligo_conc_nM Conditions supplied to [calc_tm()].
#'
#' @return An object of class `rflp_arms_primers` with `top`, `best_set`,
#'   valid-set count, parameters and exclusion diagnostics.
#' @export
design_arms_primers <- function(gene_seq, snp_pos, ref_allele, alt_allele,
                                outer_flank = 160L,
                                outer_length_min = 18L, outer_length_max = 24L,
                                inner_length_min = 18L, inner_length_max = 24L,
                                mismatch_positions = c(2L, 3L),
                                tm_min = 50, tm_max = 65, gc_min = 35, gc_max = 65,
                                dimer_dg_min = -6, hairpin_dg_min = -6,
                                heterodimer_dg_min = -6,
                                control_amplicon_min = 150L, control_amplicon_max = 500L,
                                allele_amplicon_min = 80L, allele_amplicon_max = 400L,
                                min_band_diff = 25L, max_candidates_per_pool = 20L,
                                max_raw_candidates_per_pool = Inf,
                                n_top = 5L, Na = 100, Mg = 2, dNTPs = 0.2,
                                oligo_conc_nM = 500) {
  if (max_candidates_per_pool < 1L || n_top < 1L || min_band_diff < 0L ||
      tm_min > tm_max || gc_min > gc_max || control_amplicon_min > control_amplicon_max ||
      allele_amplicon_min > allele_amplicon_max) stop("Invalid stage-2 filtering parameters.")
  pools_result <- .generate_arms_candidate_pools(
    gene_seq, snp_pos, ref_allele, alt_allele, outer_flank, outer_length_min,
    outer_length_max, inner_length_min, inner_length_max, mismatch_positions
  )
  input <- .arms_validate_inputs(gene_seq, snp_pos, ref_allele, alt_allele)
  evaluated <- lapply(pools_result$candidate_pools, .arms_evaluate_pool,
                      tm_min = tm_min, tm_max = tm_max, gc_min = gc_min, gc_max = gc_max,
                      dimer_dg_min = dimer_dg_min, hairpin_dg_min = hairpin_dg_min,
                      Na = Na, Mg = Mg, dNTPs = dNTPs, oligo_conc_nM = oligo_conc_nM,
                      max_raw_candidates = max_raw_candidates_per_pool)
  filtered_pools <- lapply(evaluated, `[[`, "candidates")
  tm_target <- mean(c(tm_min, tm_max)); gc_target <- mean(c(gc_min, gc_max))
  selected_pools <- lapply(filtered_pools, .arms_preselect_pool,
                           n_max = max_candidates_per_pool, tm_target = tm_target, gc_target = gc_target)
  diagnostics <- do.call(rbind, lapply(names(evaluated), function(name) {
    reasons <- evaluated[[name]]$reasons
    data.frame(pool = name, n_geometric = nrow(pools_result$candidate_pools[[name]]),
               n_individual_valid = nrow(filtered_pools[[name]]), n_preselected = nrow(selected_pools[[name]]),
               failed_tm = reasons[["tm"]], failed_gc = reasons[["gc"]],
               failed_self_dimer = reasons[["self_dimer"]], failed_hairpin = reasons[["hairpin"]],
               stringsAsFactors = FALSE)
  }))
  failures <- c(inner_specificity = 0L, product_size = 0L, band_separation = 0L, cross_dimer = 0L)
  valid_sets <- vector("list", 0L)
  for (layout_i in seq_len(nrow(pools_result$layouts))) {
    layout <- pools_result$layouts[layout_i, ]
    ref_pool <- selected_pools[[layout$ref_pool]]; alt_pool <- selected_pools[[layout$alt_pool]]
    if (!nrow(selected_pools$outer_forward) || !nrow(selected_pools$outer_reverse) || !nrow(ref_pool) || !nrow(alt_pool)) next
    for (f in seq_len(nrow(selected_pools$outer_forward))) for (r in seq_len(nrow(selected_pools$outer_reverse)))
      for (ri in seq_len(nrow(ref_pool))) for (ai in seq_len(nrow(alt_pool))) {
        outer_f <- selected_pools$outer_forward[f, ]; outer_r <- selected_pools$outer_reverse[r, ]
        ref_i <- ref_pool[ri, ]; alt_i <- alt_pool[ai, ]
        if (!.arms_inner_is_specific(ref_i, gene_seq, input$snp_pos, input$ref_allele, input$alt_allele) ||
            !.arms_inner_is_specific(alt_i, gene_seq, input$snp_pos, input$alt_allele, input$ref_allele)) {
          failures[["inner_specificity"]] <- failures[["inner_specificity"]] + 1L; next
        }
        control_bp <- outer_r$genomic_end - outer_f$genomic_start + 1L
        amplicon_size <- function(inner) if (identical(inner$strand, "forward"))
          outer_r$genomic_end - inner$genomic_start + 1L else inner$genomic_end - outer_f$genomic_start + 1L
        ref_bp <- amplicon_size(ref_i); alt_bp <- amplicon_size(alt_i)
        if (control_bp < control_amplicon_min || control_bp > control_amplicon_max ||
            ref_bp < allele_amplicon_min || ref_bp > allele_amplicon_max ||
            alt_bp < allele_amplicon_min || alt_bp > allele_amplicon_max) {
          failures[["product_size"]] <- failures[["product_size"]] + 1L; next
        }
        if (min(abs(c(control_bp - ref_bp, control_bp - alt_bp, ref_bp - alt_bp))) < min_band_diff) {
          failures[["band_separation"]] <- failures[["band_separation"]] + 1L; next
        }
        primers <- list(outer_f, outer_r, ref_i, alt_i)
        cross_dg <- if (is.infinite(heterodimer_dg_min) && heterodimer_dg_min < 0) rep(0, 6L) else
          vapply(utils::combn(seq_along(primers), 2L, simplify = FALSE), function(pair) {
            evaluate_dimer(primers[[pair[1]]]$sequence, primers[[pair[2]]]$sequence)$dg_approx
          }, numeric(1))
        if (any(cross_dg < heterodimer_dg_min)) { failures[["cross_dimer"]] <- failures[["cross_dimer"]] + 1L; next }
        primer_score <- sum(vapply(primers, function(p) abs(p$tm - tm_target) + abs(p$gc - gc_target) * 0.2,
                                   numeric(1)))
        row <- data.frame(layout = layout$layout,
          outer_forward_seq = outer_f$sequence, outer_forward_start = outer_f$genomic_start, outer_forward_end = outer_f$genomic_end, outer_forward_tm = outer_f$tm, outer_forward_gc = outer_f$gc,
          outer_reverse_seq = outer_r$sequence, outer_reverse_start = outer_r$genomic_start, outer_reverse_end = outer_r$genomic_end, outer_reverse_tm = outer_r$tm, outer_reverse_gc = outer_r$gc,
          ref_inner_seq = ref_i$sequence, ref_inner_start = ref_i$genomic_start, ref_inner_end = ref_i$genomic_end, ref_inner_strand = ref_i$strand, ref_inner_tm = ref_i$tm, ref_inner_gc = ref_i$gc, ref_mismatch_3prime = ref_i$mismatch_3prime,
          alt_inner_seq = alt_i$sequence, alt_inner_start = alt_i$genomic_start, alt_inner_end = alt_i$genomic_end, alt_inner_strand = alt_i$strand, alt_inner_tm = alt_i$tm, alt_inner_gc = alt_i$gc, alt_mismatch_3prime = alt_i$mismatch_3prime,
          control_amplicon_bp = control_bp, ref_amplicon_bp = ref_bp, alt_amplicon_bp = alt_bp,
          worst_cross_dimer_dg = min(cross_dg), score = round(primer_score + pmax(0, -min(cross_dg) - 4), 3), stringsAsFactors = FALSE)
        valid_sets[[length(valid_sets) + 1L]] <- row
      }
  }
  all_sets <- if (length(valid_sets)) do.call(rbind, valid_sets) else data.frame()
  if (nrow(all_sets)) {
    all_sets <- all_sets[order(all_sets$score, all_sets$layout, all_sets$outer_forward_seq, all_sets$outer_reverse_seq), , drop = FALSE]
    top <- utils::head(all_sets, n_top); top$recommended <- c("YES", rep("NO", nrow(top) - 1L)); rownames(top) <- NULL
  } else top <- data.frame()
  result <- list(top = top, best_set = if (nrow(top)) top[1, , drop = FALSE] else NULL,
                 n_geometric_sets = pools_result$n_geometric_sets, n_valid_sets = nrow(all_sets),
                 candidate_pools = pools_result$candidate_pools, filtered_pools = filtered_pools,
                 diagnostics = diagnostics, exclusion_diagnostics = as.data.frame(as.list(failures)),
                 parameters = c(pools_result$parameters, list(tm_min = tm_min, tm_max = tm_max, gc_min = gc_min, gc_max = gc_max,
                   dimer_dg_min = dimer_dg_min, hairpin_dg_min = hairpin_dg_min, heterodimer_dg_min = heterodimer_dg_min,
                   control_amplicon_min = control_amplicon_min, control_amplicon_max = control_amplicon_max,
                   allele_amplicon_min = allele_amplicon_min, allele_amplicon_max = allele_amplicon_max,
                   min_band_diff = min_band_diff, max_candidates_per_pool = max_candidates_per_pool,
                   max_raw_candidates_per_pool = max_raw_candidates_per_pool, n_top = n_top)))
  class(result) <- "rflp_arms_primers"
  result
}
