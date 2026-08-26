test_that("define_work_region retains inclusive reference coordinates", {
  reference <- arms_fixture_forward()
  region <- suppressMessages(define_work_region(
    reference,
    snp_pos = 13L,
    upstream = 6L,
    downstream = 4L
  ))

  expect_equal(region$start_region, 7L)
  expect_equal(region$end_region, 17L)
  expect_equal(region$snp_pos, 13L)
  expect_equal(length(region$sequence_region), 11L)
  expect_identical(as.character(region$sequence_region), "TTAACCGGTTA")
  expect_equal(13L - region$start_region + 1L, 7L)
})

test_that("simulate_pcr preserves coordinates and reverse-primer orientation", {
  reference <- arms_fixture_unique_amplicon()
  result <- suppressMessages(simulate_pcr(
    reference,
    fwd_primer = "TGCGTAC",
    rev_primer = "TGCAAC",
    max_mismatch = 0L
  ))

  expect_s3_class(result, "rflp_amplicon")
  expect_s4_class(result$amplicon, "DNAString")
  expect_equal(result$start, 2L)
  expect_equal(result$end, 24L)
  expect_equal(result$size, 23L)
  expect_equal(result$n_sites_fwd, 1L)
  expect_equal(result$n_sites_rev, 1L)
  expect_identical(as.character(result$amplicon), "TGCGTACCTAGGCTAACGTTGCA")
})

test_that("ARMS products retain coordinates, diagnostic sizes and allele sequence", {
  design <- design_arms_primers(
    arms_fixture_forward(), 13L, "G", "A", outer_flank = 20L,
    outer_length_min = 4L, outer_length_max = 4L, inner_length_min = 4L, inner_length_max = 4L,
    tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 100,
    dimer_dg_min = -Inf, hairpin_dg_min = -Inf, heterodimer_dg_min = -Inf,
    control_amplicon_min = 1L, control_amplicon_max = 50L,
    allele_amplicon_min = 1L, allele_amplicon_max = 50L,
    min_band_diff = 0L, max_candidates_per_pool = 3L
  )
  pcr <- simulate_arms_pcr(arms_fixture_forward(), design, 13L, "G", "A")
  heterozygote <- pcr$products[["ref/alt"]]

  expect_equal(heterozygote$size_bp, c(design$best_set$control_amplicon_bp,
                                        design$best_set$ref_amplicon_bp,
                                        design$best_set$alt_amplicon_bp))
  expect_equal(heterozygote$end - heterozygote$start + 1L, heterozygote$size_bp)
  expect_equal(heterozygote$template_allele, c("ref/alt", "ref", "alt"))
  expect_equal(substr(heterozygote$sequence[3], 13L - heterozygote$start[3] + 1L,
                      13L - heterozygote$start[3] + 1L), "A")
})
