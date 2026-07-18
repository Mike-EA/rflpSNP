# Returns TRUE if two bases are Watson-Crick complementary.
#' @keywords internal
.is_complementary <- function(b1, b2) {
  pairs <- c(A = "T", T = "A", C = "G", G = "C")
  !is.na(pairs[b1]) && pairs[b1] == b2
}

#' Evaluate dimer formation between two sequences
#'
#' Complementarity heuristic to estimate the risk of self-dimer
#' (`seq1 == seq2`) or heterodimer (`seq1 != seq2`) formation between two
#' oligonucleotides. It slides one sequence over the other at every
#' possible offset and looks for the longest run of consecutive
#' complementary bases, penalizing especially when that run occurs at the
#' 3' end (where a dimer is most critical, since it can prime polymerase
#' extension).
#'
#' @param seq1,seq2 Oligonucleotide sequences (5' -> 3'), as text or
#'   `DNAString`.
#'
#' @return A list with:
#'   \item{max_run}{Length of the longest complementary run found.}
#'   \item{dg_approx}{Heuristic estimate of free energy (kcal/mol; more
#'     negative values indicate higher risk). **Not a rigorous
#'     thermodynamic calculation**; useful only for comparing candidates
#'     against each other.}
#'   \item{risk_3prime}{`TRUE` if the longest run involves the 3' end of
#'     `seq1`.}
#'
#' @details
#' This is a fast screening heuristic, not a substitute for a rigorous
#' thermodynamic calculator (e.g. IDT OligoAnalyzer or Primer3). Verifying
#' the final pair before synthesis is recommended.
#'
#' @export
evaluate_dimer <- function(seq1, seq2) {
  s1 <- strsplit(toupper(as.character(seq1)), "")[[1]]
  s2 <- strsplit(toupper(as.character(seq2)), "")[[1]]
  s2_rev <- rev(s2)
  n1 <- length(s1)
  n2 <- length(s2_rev)

  best_run <- 0
  run_at_3prime <- FALSE

  for (offset in -(n2 - 1):(n1 - 1)) {
    current_run <- 0
    local_max_run <- 0
    final_run_pos <- NA
    for (i in seq_len(n1)) {
      j <- i - offset
      if (j >= 1 && j <= n2 && .is_complementary(s1[i], s2_rev[j])) {
        current_run <- current_run + 1
        if (current_run > local_max_run) {
          local_max_run <- current_run
          final_run_pos <- i
        }
      } else {
        current_run <- 0
      }
    }
    if (local_max_run > best_run) {
      best_run <- local_max_run
      run_at_3prime <- !is.na(final_run_pos) && final_run_pos >= (n1 - 4)
    }
  }

  dg_approx <- -1.5 * best_run
  if (run_at_3prime) dg_approx <- dg_approx * 1.3

  list(max_run = best_run, dg_approx = round(dg_approx, 1), risk_3prime = run_at_3prime)
}

#' Evaluate hairpin (intramolecular) structure formation
#'
#' Heuristic to detect hairpin (stem-loop) structures within a single
#' primer: it searches for the longest pair of internal arms that are
#' complementary to each other (one read 5'->3' and the other 3'->5'),
#' separated by a loop of at least `min_loop` nucleotides. This structure
#' is not evaluated by [evaluate_dimer()], which only compares two
#' sequences against each other (intermolecular), not a sequence against
#' itself in a hairpin configuration (intramolecular).
#'
#' @param primer_seq Primer sequence (5' -> 3'), as text or `DNAString`.
#' @param min_stem Minimum stem length to consider a risk (default `3`).
#' @param min_loop Minimum loop length between the two stem arms (default
#'   `3`; a real hairpin requires a minimum loop size to fold physically).
#'
#' @return A list with:
#'   \item{stem}{Length of the longest complementary stem found.}
#'   \item{loop}{Length of the loop associated with that stem (`NA` if
#'     `stem == 0`).}
#'   \item{dg_approx}{Heuristic estimate of free energy (kcal/mol; see the
#'     note in [evaluate_dimer()] about its limitations).}
#'   \item{risk}{`TRUE` if `stem >= 4`, the default heuristic threshold for
#'     flagging the hairpin as concerning.}
#'
#' @export
evaluate_hairpin <- function(primer_seq, min_stem = 3, min_loop = 3) {
  pair_map <- c(A = "T", T = "A", C = "G", G = "C")
  s <- strsplit(toupper(as.character(primer_seq)), "")[[1]]
  n <- length(s)

  best_stem <- 0
  best_loop <- NA_integer_

  if (n < (2 * min_stem + min_loop)) {
    return(list(stem = 0, loop = NA_integer_, dg_approx = 0, risk = FALSE))
  }

  max_possible_loop <- n - 2 * min_stem
  for (loop_size in min_loop:max_possible_loop) {
    max_arm <- floor((n - loop_size) / 2)
    if (max_arm < min_stem) next

    for (arm_len in min_stem:max_arm) {
      max_start5 <- n - loop_size - 2 * arm_len + 1
      if (max_start5 < 1) next

      for (start5 in 1:max_start5) {
        arm5 <- s[start5:(start5 + arm_len - 1)]
        start3 <- start5 + arm_len + loop_size
        end3 <- start3 + arm_len - 1
        if (end3 > n) next

        arm3 <- s[start3:end3]
        arm3_rev <- rev(arm3)

        complementary_count <- sum(mapply(
          function(b1, b2) .is_complementary(b1, b2), arm5, arm3_rev
        ))

        if (complementary_count == arm_len && arm_len > best_stem) {
          best_stem <- arm_len
          best_loop <- loop_size
        }
      }
    }
  }

  dg_approx <- -1.8 * best_stem
  risk <- best_stem >= 4

  list(stem = best_stem, loop = best_loop, dg_approx = round(dg_approx, 1), risk = risk)
}
