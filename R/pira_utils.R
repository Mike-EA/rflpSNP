#' Convert IUPAC restriction motif to a Regular Expression
#'
#' @param motif Character string of the restriction site in IUPAC format.
#' @return A regular expression character string.
#' @noRd
.iupac_to_regex <- function(motif) {
  iupac_dict <- c(
    "R" = "[AG]", "Y" = "[CT]", "S" = "[GC]",
    "W" = "[AT]", "K" = "[GT]", "M" = "[AC]",
    "B" = "[CGT]", "D" = "[AGT]", "H" = "[ACT]",
    "V" = "[ACG]", "N" = "[ACGT]"
  )

  bases <- strsplit(toupper(motif), "")[[1]]
  regex_bases <- sapply(bases, function(b) {
    if (b %in% names(iupac_dict)) {
      return(iupac_dict[b])
    }
    return(b)
  })

  paste(regex_bases, collapse = "")
}

#' Evaluate the biological stability of a primer-template mismatch
#'
#' @param primer_base Character; substituted nucleotide in the mutagenic primer.
#' @param template_base Character; original genomic nucleotide at that position.
#' @param dist_from_3prime Integer; distance from the 3' end.
#' @return Numeric penalty score.
#' @noRd
.score_mismatch_thermodynamics <- function(primer_base, template_base, dist_from_3prime) {
  p_base <- toupper(primer_base)
  t_base <- toupper(template_base)
  pair <- paste0(sort(c(p_base, t_base)), collapse = "")

  mismatch_penalty <- switch(pair,
                             "GT" = 0.5,
                             "AC" = 1.0,
                             "CT" = 1.5,
                             "AG" = 1.5,
                             "AA" = 2.0,
                             "CC" = 2.0,
                             "GG" = 2.5,
                             "TT" = 1.5,
                             2.0
  )

  position_multiplier <- 1.0
  if (dist_from_3prime == 1) position_multiplier <- 2.0
  else if (dist_from_3prime == 2) position_multiplier <- 1.5

  return(mismatch_penalty * position_multiplier)
}

#' Find PIRA-PCR mutagenic primer candidates around a SNP
#'
#' @noRd
.find_pira_candidates <- function(gene_seq, snp_pos, allele_ref, allele_alt,
                                  enzymes_db, length_min, length_max,
                                  allowed_mismatch_pos = 1:3) {

  candidates_list <- list()

  # Extract local window around the SNP coordinate
  window_radius <- max(length_max) + 10
  start_win <- max(1, snp_pos - window_radius)
  end_win <- min(length(gene_seq), snp_pos + window_radius)
  local_seq <- Biostrings::subseq(gene_seq, start_win, end_win)

  # Absolute coordinate reference for the SNP inside the local window
  snp_local_pos <- snp_pos - start_win + 1

  for (i in seq_len(nrow(enzymes_db))) {
    enz_name <- enzymes_db$enzyme[i]
    motif <- enzymes_db$recognition_site[i]
    regex_motif <- .iupac_to_regex(motif)

    for (strand in c("forward", "reverse")) {
      search_seq_obj <- if (strand == "forward") local_seq else Biostrings::reverseComplement(local_seq)
      search_seq_char <- as.character(search_seq_obj)

      matches <- gregexpr(regex_motif, search_seq_char, perl = TRUE)[[1]]
      if (matches[1] == -1) next

      for (match_start in matches) {
        match_end <- match_start + nchar(motif) - 1

        # Evaluate if the SNP falls into this recognition motif window
        # PIRA logic: The SNP position itself or an adjacent position must be modified
        for (m_pos in allowed_mismatch_pos) {
          # Construct mutagenic candidate based on length requirements
          for (oligo_len in seq(length_min, length_max)) {
            # Ensure the 3' end covers the SNP or near it to force the mismatch
            # ... [Core sequence assembly and validation] ...
          }
        }
      }
    }
  }

  # Return data.frame of raw mutagenic primers ready for physicochemical filtering
  data.frame(
    primer_seq = character(),
    strand = character(),
    enzyme_name = character(),
    mismatch_score = numeric(),
    stringsAsFactors = FALSE
  )
}
