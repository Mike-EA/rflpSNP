test_that("FTO rs8050136 regression case returns documented products, report and gel", {
  reference <- arms_fixture_fto_rs8050136()
  snp <- suppressMessages(locate_snp(reference, "CATGCCAGTTGCCCACTGTGGCAAT"))
  expect_equal(snp$snp_pos, 86L)
  expect_equal(snp$snp_base, "C")

  design <- design_arms_primers(
    reference, snp$snp_pos, "C", "A",
    outer_flank = 85L, outer_length_min = 20L, outer_length_max = 20L,
    inner_length_min = 18L, inner_length_max = 20L,
    tm_min = 55, tm_max = 72, gc_min = 35, gc_max = 65,
    dimer_dg_min = -8, hairpin_dg_min = -6, heterodimer_dg_min = -9,
    control_amplicon_min = 150L, control_amplicon_max = 200L,
    allele_amplicon_min = 80L, allele_amplicon_max = 150L,
    min_band_diff = 20L, max_candidates_per_pool = 8L,
    max_raw_candidates_per_pool = 8L
  )
  # Candidate counts may vary slightly with numerical dependency versions;
  # the recommended set and its simulated products are the regression contract.
  expect_true(design$n_valid_sets >= 1L)
  expect_equal(unname(unlist(design$best_set[c("control_amplicon_bp", "ref_amplicon_bp", "alt_amplicon_bp")])), c(152, 105, 84))
  expect_equal(design$best_set$outer_forward_seq, "GTATTTGATTTCCTTTTCCC")
  expect_equal(design$best_set$outer_reverse_seq, "TTTCCAAGCATTCCATGAGT")

  pcr <- simulate_arms_pcr(reference, design, snp$snp_pos, "C", "A")
  expect_equal(pcr$products[["ref/alt"]]$size_bp, c(152, 105, 84))
  report <- tempfile(fileext = ".txt")
  expect_message(export_arms_primers_txt(design, report), "exported")
  expect_match(paste(readLines(report, warn = FALSE), collapse = "\n"), "152 bp")
  expect_s3_class(simulate_arms_gel(pcr), "ggplot")

  map_file <- tempfile(fileext = ".txt")
  expect_message(export_arms_amplicon_map_txt(pcr, map_file, line_width = 40L), "map exported")
  map_text <- paste(readLines(map_file, warn = FALSE), collapse = "\n")
  expect_match(map_text, "GENOTYPE: ref/alt")
  expect_match(map_text, "GTATTTGATTTCCTTTTCCC")
  expect_match(map_text, "Complement3'")
})
