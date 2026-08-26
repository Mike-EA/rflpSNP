test_that("TmCalculator backend returns a valid Tm", {
  expect_gte(
    utils::packageVersion("TmCalculator"),
    numeric_version("1.0.8")
  )

  tm <- calc_tm("ATGCGATGCGATGCATGCA")

  expect_type(tm, "double")
  expect_length(tm, 1)
  expect_true(is.finite(tm))
  expect_equal(tm, 70.5, tolerance = 0.2)
})

test_that("check_tm_backend reports success and returns the test value", {
  output <- utils::capture.output(
    result <- check_tm_backend()
  )

  expect_match(paste(output, collapse = "\n"), "\\[OK\\]")
  expect_true(is.numeric(result))
  expect_true(is.finite(result))
})
