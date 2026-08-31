test_that("ARMS synthetic fixtures are DNAString objects with fixed coordinates", {
  expect_s4_class(arms_fixture_forward(), "DNAString")
  expect_equal(length(arms_fixture_forward()), 20)
  expect_equal(as.character(arms_fixture_forward()[13]), "G")

  expect_s4_class(arms_fixture_reverse(), "DNAString")
  expect_equal(length(arms_fixture_reverse()), 19)
  expect_equal(as.character(arms_fixture_reverse()[5]), "A")
})

test_that("locate_snp uses 1-based coordinates on the forward reference", {
  result <- suppressMessages(locate_snp(
    arms_fixture_forward(),
    flank_seq = "GGTTAACC"
  ))

  expect_equal(result$snp_pos, 13L)
  expect_equal(result$snp_base, "G")
  expect_identical(result$strand, "forward")
})

test_that("locate_snp translates a complementary-strand flank to the reference", {
  result <- suppressMessages(locate_snp(
    arms_fixture_reverse(),
    flank_seq = "GATCGT"
  ))

  expect_equal(result$snp_pos, 5L)
  expect_equal(result$snp_base, "A")
  expect_identical(result$strand, "complementary")
})

test_that("ARMS conventions retain reference coordinates and 3-prime alleles", {
  reference <- arms_fixture_forward()
  snp_pos <- 13L
  ref_allele <- as.character(reference[snp_pos])
  alt_allele <- "A"

  expect_identical(ref_allele, "G")
  expect_match(ref_allele, "^[ACGT]$")
  expect_match(alt_allele, "^[ACGT]$")
  expect_false(identical(ref_allele, alt_allele))

  forward_internal <- Biostrings::subseq(reference, 7L, snp_pos)
  reverse_internal <- Biostrings::reverseComplement(
    Biostrings::subseq(reference, snp_pos, 19L)
  )

  expect_identical(as.character(forward_internal[length(forward_internal)]), ref_allele)
  expect_identical(
    as.character(reverse_internal[length(reverse_internal)]),
    as.character(Biostrings::complement(Biostrings::DNAString(ref_allele)))
  )
  expect_identical(
    as.character(Biostrings::reverseComplement(reverse_internal)),
    as.character(Biostrings::subseq(reference, snp_pos, 19L))
  )
  expect_equal(length(forward_internal), snp_pos - 7L + 1L)
  expect_equal(length(reverse_internal), 19L - snp_pos + 1L)
})

test_that("locate_snp documents ambiguous and invalid synthetic inputs", {
  expect_warning(
    repeated <- suppressMessages(locate_snp(
      arms_fixture_repeated_flank(),
      flank_seq = "AACCGG"
    )),
    "matches 2 times"
  )
  expect_equal(repeated$snp_pos, 7L)

  expect_error(
    locate_snp(Biostrings::DNAString("AACCGG"), "AACCGG"),
    "outside the sequence bounds"
  )
  expect_error(
    locate_snp(arms_fixture_forward(), "TTTTTT"),
    "not found on either strand"
  )
})
