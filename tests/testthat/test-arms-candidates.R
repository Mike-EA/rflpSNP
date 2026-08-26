test_that("ARMS candidate pools enforce geometry, alleles and deliberate mismatches", {
  candidates <- rflpSNP:::.generate_arms_candidate_pools(
    arms_fixture_forward(),
    snp_pos = 13L,
    ref_allele = "G",
    alt_allele = "A",
    outer_flank = 20L,
    outer_length_min = 4L,
    outer_length_max = 4L,
    inner_length_min = 4L,
    inner_length_max = 4L
  )

  expect_s3_class(candidates, "rflp_arms_primers")
  expect_true(is.data.frame(candidates$top))
  expect_null(candidates$best_set)
  expect_equal(nrow(candidates$candidate_pools$outer_forward), 9L)
  expect_equal(nrow(candidates$candidate_pools$outer_reverse), 4L)
  expect_true(all(candidates$candidate_pools$outer_forward$genomic_end < 13L))
  expect_true(all(candidates$candidate_pools$outer_reverse$genomic_start > 13L))

  expect_true(all(candidates$candidate_pools$ref_forward$genomic_end == 13L))
  expect_true(all(candidates$candidate_pools$ref_reverse$genomic_start == 13L))
  expect_true(all(substr(candidates$candidate_pools$ref_forward$sequence, 4L, 4L) == "G"))
  expect_true(all(substr(candidates$candidate_pools$alt_forward$sequence, 4L, 4L) == "A"))
  expect_true(all(substr(candidates$candidate_pools$ref_reverse$sequence, 4L, 4L) == "C"))
  expect_true(all(substr(candidates$candidate_pools$alt_reverse$sequence, 4L, 4L) == "T"))
  expect_setequal(candidates$candidate_pools$ref_forward$mismatch_3prime, c("3'-2", "3'-3"))
  expect_true(all(candidates$candidate_pools$ref_forward$mismatch_from != candidates$candidate_pools$ref_forward$mismatch_to))

  expect_equal(candidates$n_geometric_sets, 2592)
  expect_setequal(candidates$layouts$layout, c("ref_forward_alt_reverse", "alt_forward_ref_reverse"))
  expect_equal(sum(candidates$diagnostics$n_candidates), 37L)
  expect_match(paste(capture.output(print(candidates)), collapse = "\n"), "Filtering and ranking")
})

test_that("ARMS candidate pools reject an inconsistent reference allele", {
  expect_error(
    rflpSNP:::.generate_arms_candidate_pools(
      arms_fixture_forward(), 13L, ref_allele = "A", alt_allele = "G",
      outer_length_min = 4L, outer_length_max = 4L,
      inner_length_min = 4L, inner_length_max = 4L
    ),
    "does not match"
  )
})

test_that("ARMS candidate pools report impossible geometry and reject bad mismatch positions", {
  impossible <- rflpSNP:::.generate_arms_candidate_pools(
    arms_fixture_forward(), 2L, ref_allele = "A", alt_allele = "C",
    outer_flank = 1L, outer_length_min = 8L, outer_length_max = 8L,
    inner_length_min = 8L, inner_length_max = 8L
  )
  expect_equal(impossible$n_geometric_sets, 0)
  expect_equal(impossible$diagnostics$n_candidates[impossible$diagnostics$pool == "outer_forward"], 0L)
  expect_equal(impossible$diagnostics$n_candidates[impossible$diagnostics$pool == "outer_reverse"], 0L)
  expect_match(paste(capture.output(print(impossible)), collapse = "\n"), "No complete")

  expect_error(
    rflpSNP:::.generate_arms_candidate_pools(
      arms_fixture_forward(), 13L, ref_allele = "G", alt_allele = "A",
      outer_length_min = 4L, outer_length_max = 4L,
      inner_length_min = 4L, inner_length_max = 4L,
      mismatch_positions = 1L
    ),
    "Invalid candidate"
  )
})

test_that("ARMS geometry count remains numeric for large candidate pools", {
  reference <- Biostrings::DNAString(paste(rep("ACGT", 1250L), collapse = ""))
  candidates <- rflpSNP:::.generate_arms_candidate_pools(
    reference, 2501L, ref_allele = "A", alt_allele = "C",
    outer_flank = 2000L, outer_length_min = 1L, outer_length_max = 1L,
    inner_length_min = 18L, inner_length_max = 80L, mismatch_positions = 2L
  )
  expect_type(candidates$n_geometric_sets, "double")
  expect_gt(candidates$n_geometric_sets, .Machine$integer.max)
  expect_false(is.na(candidates$n_geometric_sets))
})

test_that("design_arms_primers filters, ranks and reports exclusions reproducibly", {
  design <- design_arms_primers(
    arms_fixture_forward(), snp_pos = 13L, ref_allele = "G", alt_allele = "A",
    outer_flank = 20L, outer_length_min = 4L, outer_length_max = 4L,
    inner_length_min = 4L, inner_length_max = 4L,
    tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 100,
    dimer_dg_min = -Inf, hairpin_dg_min = -Inf, heterodimer_dg_min = -Inf,
    control_amplicon_min = 1L, control_amplicon_max = 50L,
    allele_amplicon_min = 1L, allele_amplicon_max = 50L,
    min_band_diff = 0L, max_candidates_per_pool = 3L, n_top = 5L
  )

  expect_s3_class(design, "rflp_arms_primers")
  expect_gt(design$n_valid_sets, 0L)
  expect_equal(nrow(design$top), 5L)
  expect_equal(design$top$recommended[1], "YES")
  expect_true(is.data.frame(design$best_set))
  expect_true(all(c("failed_tm", "failed_gc", "failed_self_dimer", "failed_hairpin") %in% names(design$diagnostics)))
  expect_true(all(c("inner_specificity", "product_size", "band_separation", "cross_dimer") %in% names(design$exclusion_diagnostics)))
  expect_true(all(diff(design$top$score) >= 0))
})

test_that("design_arms_primers retains diagnostics when no set survives", {
  design <- design_arms_primers(
    arms_fixture_forward(), snp_pos = 13L, ref_allele = "G", alt_allele = "A",
    outer_flank = 20L, outer_length_min = 4L, outer_length_max = 4L,
    inner_length_min = 4L, inner_length_max = 4L,
    tm_min = 1000, tm_max = 1001
  )
  expect_equal(design$n_valid_sets, 0L)
  expect_equal(nrow(design$top), 0L)
  expect_null(design$best_set)
  expect_gt(sum(design$diagnostics$failed_tm), 0L)
  expect_match(paste(capture.output(print(design)), collapse = "\\n"), "No set passed")
})
