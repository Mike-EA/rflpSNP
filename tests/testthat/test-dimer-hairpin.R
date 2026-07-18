test_that("evaluate_dimer detects a strong self-dimer", {
  r <- evaluate_dimer("GGGGGCCCCC", "GGGGGCCCCC")
  expect_equal(r$max_run, 10)
  expect_true(r$dg_approx < -15)
  expect_true(r$risk_3prime)
})

test_that("evaluate_dimer reports no risk when there is no complementarity", {
  r <- evaluate_dimer("AAAAAAAAAA", "AAAAAAAAAA")
  expect_equal(r$max_run, 0)
  expect_equal(r$dg_approx, 0)
  expect_false(r$risk_3prime)
})

test_that("evaluate_hairpin detects a clear hairpin (stem=5, loop=3)", {
  h <- evaluate_hairpin("GGGGGAAACCCCC")
  expect_equal(h$stem, 5)
  expect_equal(h$loop, 3)
  expect_true(h$risk)
})

test_that("evaluate_hairpin does not fail on very short sequences", {
  h <- evaluate_hairpin("ATCG")
  expect_equal(h$stem, 0)
  expect_false(h$risk)
})

test_that("evaluate_dimer and evaluate_hairpin assess different things", {
  # A primer with no relevant self-dimer (sequence against itself) may
  # still carry hairpin risk if it contains an internal inverted stem.
  primer <- "ATGCGGGGGAAACCCCCATGC"
  dimer <- evaluate_dimer(primer, primer)
  hairpin <- evaluate_hairpin(primer)
  expect_true(hairpin$stem >= 5)
  expect_true(hairpin$risk)
})
