# --- Hand-verified synthetic case -------------------------------------
# Built so that EcoRI (GAATTC) can only be completed, for the REFERENCE
# allele, by introducing exactly 1 artificial mismatch at position
# snp_pos - 1 (forward orientation). Worked out by hand in the Fase 3
# planning notes; see comments below for the position-by-position trace.
#
# snp_pos = 10. EcoRI site spans positions 6-11 when snp_pos sits at
# offset 5 (the first "T" in G-A-A-T-T-C):
#   pos  6  7  8  9 10 11   <- genomic coordinate
#   code G  A  A  T [T] C   <- EcoRI pattern (snp_pos in brackets)
#   seq  G  A  A  G  T  C   <- what we actually put in gene_seq
# Position 9 (code T) has genomic "G" -> needs 1 artificial mismatch.
# Positions 6, 7, 8, 11 already match the pattern with no changes needed.
# reference_allele = "T" completes the site; alternate_allele = "C" does not.

build_synthetic_seq <- function() {
  Biostrings::DNAString("ACGTAGAAGTCACGTA")
}

test_that("scan_pira_candidates finds the expected hand-verified EcoRI candidate", {
  gene_seq <- build_synthetic_seq()
  snp_pos <- 10L

  result <- scan_pira_candidates(
    gene_seq, snp_pos,
    reference_allele = "T", alternate_allele = "C",
    enzymes = restriction_enzymes[restriction_enzymes$enzyme == "EcoRI", ],
    orientations = "forward"
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$enzyme, "EcoRI")
  expect_equal(result$orientation, "forward")
  expect_equal(result$site_start, 6L)
  expect_equal(result$site_end, 11L)
  expect_equal(result$n_artificial_mismatches, 1L)
  expect_equal(result$mismatch_positions, "9")
  expect_equal(result$mismatch_bases_required, "T")
  expect_equal(result$creates_site_for, "reference")
})

test_that("scan_pira_candidates returns nothing when alleles both fail to complete the site", {
  gene_seq <- build_synthetic_seq()
  snp_pos <- 10L

  # Neither "A" nor "G" is "T" -> neither completes the EcoRI site at
  # this position; not reported (no discrimination possible).
  result <- scan_pira_candidates(
    gene_seq, snp_pos,
    reference_allele = "A", alternate_allele = "G",
    enzymes = restriction_enzymes[restriction_enzymes$enzyme == "EcoRI", ],
    orientations = "forward"
  )
  expect_equal(nrow(result), 0)
})

test_that("scan_pira_candidates returns a well-formed empty result when nothing is achievable within the mismatch budget", {
  gene_seq <- Biostrings::DNAString(paste0(rep("C", 20), collapse = ""))
  fake_enzyme <- data.frame(
    enzyme = "TestEnz7A",
    recognition_site = "AAAAAAA",  # 7 bp, all fixed 'A' -- genomic is all 'C'
    cut_type = "blunt",
    stringsAsFactors = FALSE
  )
  result <- scan_pira_candidates(
    gene_seq, snp_pos = 10L,
    reference_allele = "A", alternate_allele = "C",
    enzymes = fake_enzyme,
    orientations = c("forward", "reverse")
  )
  expect_equal(nrow(result), 0)
  expect_named(result, c("enzyme", "orientation", "site_start", "site_end",
                          "n_artificial_mismatches", "mismatch_positions",
                          "mismatch_bases_required", "creates_site_for",
                          "cut_type", "recognition_site"))
})

test_that("scan_pira_candidates validates its inputs", {
  gene_seq <- build_synthetic_seq()
  expect_error(scan_pira_candidates(gene_seq, snp_pos = 100L,
                                     reference_allele = "T", alternate_allele = "C"),
               "within the bounds")
  expect_error(scan_pira_candidates(gene_seq, snp_pos = 10L,
                                     reference_allele = "T", alternate_allele = "T"),
               "must differ")
  expect_error(scan_pira_candidates(gene_seq, snp_pos = 10L,
                                     reference_allele = "X", alternate_allele = "C"),
               "must each be one of A, C, G, T")
})

# --- Live integration test: real locus, rs8050136 (FTO, chr16:53782363) ----
# Fetches ~100 bp of real GRCh38 sequence around the SNP directly from
# Ensembl and runs the full package pipeline against it. Confirmed
# alleles (2026-07-21, NCBI dbSNP build 157): C (reference) / A (alternate),
# both on the + strand, consistent with the FTO gene's + strand orientation
# -- no strand-flip expected here (unlike rs1801133).

test_that("scan_pira_candidates runs against the real rs8050136 (FTO) locus", {
  fetch_region <- function() {
    resp <- httr::GET(
      "https://rest.ensembl.org/sequence/region/human/16:53782313-53782413:1",
      httr::accept("text/x-fasta"), httr::timeout(20)
    )
    if (httr::status_code(resp) != 200) {
      stop("Ensembl sequence request failed with status ", httr::status_code(resp))
    }
    httr::content(resp, as = "text", encoding = "UTF-8")
  }

  # Self-contained network guard (does not depend on helper-network.R being
  # loaded first -- see commit notes on an observed, unexplained
  # helper-loading order issue affecting newly-added test files).
  region <- tryCatch(
    fetch_region(),
    error = function(e) {
      msg <- conditionMessage(e)
      network_pattern <- paste("internet connection", "Timeout", "timed out",
                                "Could not reach", "Failed to connect",
                                "status", sep = "|")
      if (grepl(network_pattern, msg, ignore.case = TRUE)) {
        skip(paste("Network/API issue, skipping test:", msg))
      }
      stop(e)
    }
  )

  fasta_lines <- strsplit(region, "\n")[[1]]
  seq_text <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
  gene_seq <- Biostrings::DNAString(seq_text)

  # 53782363 is the 51st base of the fetched 53782313-53782413 window.
  snp_pos <- 53782363L - 53782313L + 1L
  expect_equal(as.character(gene_seq[snp_pos]), "C")  # sanity check vs dbSNP

  result <- scan_pira_candidates(
    gene_seq, snp_pos,
    reference_allele = "C", alternate_allele = "A"
  )

  # We don't assert a specific enzyme here (unlike the synthetic case,
  # this is real, previously-unexamined sequence) -- just that the
  # pipeline runs end-to-end on real data and returns a well-formed
  # result, and that IF candidates are found, they're internally
  # consistent (each site actually spans snp_pos, orientation/columns
  # are valid, mismatch count respects the budget).
  expect_true(is.data.frame(result))
  if (nrow(result) > 0) {
    expect_true(all(result$site_start <= snp_pos & result$site_end >= snp_pos))
    expect_true(all(result$n_artificial_mismatches %in% 1:2))
    expect_true(all(result$creates_site_for %in% c("reference", "alternate")))
    expect_true(all(result$orientation %in% c("forward", "reverse")))
  }
})
