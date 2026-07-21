# NOTE: Mike had no existing tests for locate_snp(). These cover both the
# pre-existing (unchanged) behavior and the new optional `rsid` parameter.
# The unit tests below use short synthetic sequences and require no network
# access; the integration test at the bottom requires internet and reuses
# the shared `run_or_skip_offline()` helper (tests/testthat/helper-network.R).

flank <- "GAAAAGCTGCGTGATGATGAAATCG"  # same flank used in the original MTHFR script

test_that("locate_snp finds a forward-strand match and reads the correct base", {
  gene_seq <- Biostrings::DNAString(paste0("ACGT", flank, "C", "ACGTACGT"))

  result <- locate_snp(gene_seq, flank_seq = flank)

  expect_equal(result$strand, "forward")
  expect_equal(result$snp_base, "C")
  expect_equal(result$snp_pos, 4 + nchar(flank) + 1)
  expect_null(result$alleles)  # rsid not supplied -> unchanged return shape
})

test_that("locate_snp finds a complementary-strand match and returns a self-consistent base", {
  forward_construction <- Biostrings::DNAString(paste0("ACGT", flank, "C", "ACGTACGT"))
  gene_seq <- Biostrings::reverseComplement(forward_construction)

  result <- locate_snp(gene_seq, flank_seq = flank)

  expect_equal(result$strand, "complementary")
  # Self-consistency: whatever position is reported, re-reading gene_seq at
  # that position must give exactly the returned snp_base.
  expect_equal(result$snp_base, as.character(gene_seq[result$snp_pos]))
})

test_that("locate_snp warns when the flank matches more than once", {
  gene_seq <- Biostrings::DNAString(paste0("ACGT", flank, "C", flank, "G"))
  expect_warning(
    result <- locate_snp(gene_seq, flank_seq = flank),
    "matches 2 times on the forward strand"
  )
  expect_equal(result$strand, "forward")
})

test_that("locate_snp errors when the flank is not found on either strand", {
  gene_seq <- Biostrings::DNAString("ACGTACGTACGTACGTACGTACGTACGTACGT")
  expect_error(
    locate_snp(gene_seq, flank_seq = flank),
    "not found on either strand"
  )
})

test_that("locate_snp errors when the computed SNP position falls out of bounds", {
  # Flank ends exactly at the last base of gene_seq -> default snp_offset = 1
  # points one base past the end.
  gene_seq <- Biostrings::DNAString(flank)
  expect_error(
    locate_snp(gene_seq, flank_seq = flank),
    "falls outside the sequence bounds"
  )
})

test_that("locate_snp with rsid resolves the correct PIRA-PCR alternate alleles for MTHFR C677T", {
  gene_seq <- Biostrings::DNAString(paste0("ACGT", flank, "C", "ACGTACGT"))

  # Self-contained network guard (does not depend on helper-network.R being
  # loaded first -- see note in PROVENANCE/commit message about an observed
  # helper-loading order issue).
  result <- tryCatch(
    locate_snp(gene_seq, flank_seq = flank, rsid = "rs1801133"),
    error = function(e) {
      msg <- conditionMessage(e)
      network_pattern <- paste("internet connection", "Timeout", "timed out",
                                "Could not reach", "Failed to connect", sep = "|")
      if (grepl(network_pattern, msg, ignore.case = TRUE)) {
        skip(paste("Network/API issue, skipping test:", msg))
      }
      stop(e)
    }
  )

  # Unchanged core behavior:
  expect_equal(result$strand, "forward")
  expect_equal(result$snp_base, "C")

  # New: allele resolution attached under result$alleles
  expect_false(is.null(result$alleles))
  expect_equal(result$alleles$reference_allele, "C")
  expect_equal(result$alleles$strand_used, "complement")
  expect_true(result$alleles$orientation_verified)
  expect_true(result$alleles$multiallelic)
  expect_setequal(result$alleles$alternate_alleles, c("T", "G"))
})
