test_that("locate_snp explains complementary-strand results", {
  gene_seq <- Biostrings::DNAString("TTTACGATCTAAA")

  output <- utils::capture.output(
    result <- locate_snp(gene_seq, flank_seq = "GATCGT"),
    type = "message"
  )

  expect_true(any(grepl("equivalent base in the flank orientation is A", output, fixed = TRUE)))
  expect_equal(result$snp_pos, 3)
  expect_identical(result$snp_base, "T")
  expect_identical(result$strand, "complementary")
})

test_that("locate_snp does not print the complementary-strand note for a forward match", {
  gene_seq <- Biostrings::DNAString("TTTGATCGTAAA")

  output <- utils::capture.output(
    result <- locate_snp(gene_seq, flank_seq = "GATCGT"),
    type = "message"
  )

  expect_false(any(grepl("equivalent base in the flank orientation", output, fixed = TRUE)))
  expect_equal(result$snp_pos, 10)
  expect_identical(result$snp_base, "A")
  expect_identical(result$strand, "forward")
})
