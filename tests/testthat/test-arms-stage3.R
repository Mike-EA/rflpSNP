arms_stage2_design <- function() {
  design_arms_primers(
    arms_fixture_forward(), snp_pos = 13L, ref_allele = "G", alt_allele = "A",
    outer_flank = 20L, outer_length_min = 4L, outer_length_max = 4L,
    inner_length_min = 4L, inner_length_max = 4L,
    tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 100,
    dimer_dg_min = -Inf, hairpin_dg_min = -Inf, heterodimer_dg_min = -Inf,
    control_amplicon_min = 1L, control_amplicon_max = 50L,
    allele_amplicon_min = 1L, allele_amplicon_max = 50L,
    min_band_diff = 0L, max_candidates_per_pool = 3L, n_top = 5L
  )
}

test_that("simulate_arms_pcr returns the expected bands for three genotypes", {
  design <- arms_stage2_design()
  simulation <- simulate_arms_pcr(arms_fixture_forward(), design, 13L, "G", "A")

  expect_s3_class(simulation, "rflp_arms_pcr")
  expect_equal(names(simulation$products), c("ref/ref", "ref/alt", "alt/alt"))
  expect_equal(simulation$products[["ref/ref"]]$product, c("control", "reference"))
  expect_equal(simulation$products[["ref/alt"]]$product, c("control", "reference", "alternative"))
  expect_equal(simulation$products[["alt/alt"]]$product, c("control", "alternative"))
  expect_equal(simulation$products[["ref/alt"]]$size_bp,
               c(design$best_set$control_amplicon_bp, design$best_set$ref_amplicon_bp, design$best_set$alt_amplicon_bp))
  alt_product <- simulation$products[["alt/alt"]][simulation$products[["alt/alt"]]$product == "alternative", ]
  expect_equal(substr(alt_product$sequence, 13L - alt_product$start + 1L, 13L - alt_product$start + 1L), "A")

  broken <- design$best_set
  broken$ref_inner_end <- broken$ref_inner_end - 1L
  expect_error(simulate_arms_pcr(arms_fixture_forward(), broken, 13L, "G", "A"), "Inner-primer")
})

test_that("ARMS report, gel and pipeline expose all stage-3 outputs", {
  design <- arms_stage2_design()
  report <- tempfile(fileext = ".txt")
  expect_message(export_arms_primers_txt(design, report), "exported")
  expect_true(file.exists(report))
  report_text <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(report_text, "TETRA-PRIMER ARMS-PCR")
  expect_match(report_text, "ref/alt")

  simulation <- simulate_arms_pcr(arms_fixture_forward(), design$best_set, 13L, "G", "A")
  expect_s3_class(simulate_arms_gel(simulation), "ggplot")

  pipeline <- suppressMessages(run_arms_pcr_pipeline(
    arms_fixture_forward(), flank_seq = "GGTTAACC", alt_allele = "A", export_txt = FALSE,
    outer_flank = 20L, outer_length_min = 4L, outer_length_max = 4L,
    inner_length_min = 4L, inner_length_max = 4L,
    tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 100,
    dimer_dg_min = -Inf, hairpin_dg_min = -Inf, heterodimer_dg_min = -Inf,
    control_amplicon_min = 1L, control_amplicon_max = 50L,
    allele_amplicon_min = 1L, allele_amplicon_max = 50L,
    min_band_diff = 0L, max_candidates_per_pool = 3L
  ))
  expect_named(pipeline, c("snp", "design", "report", "pcr_result", "gel"))
  expect_null(pipeline$report)
  expect_s3_class(pipeline$gel, "ggplot")
})
