test_that("fetch_dbsnp_alleles validates rsid format without network access", {
  expect_error(fetch_dbsnp_alleles(123), "must be a single character string")
  expect_error(fetch_dbsnp_alleles("not_an_rsid"), "does not look like a valid")
  expect_error(fetch_dbsnp_alleles("RS12345"), "does not look like a valid")
  expect_error(fetch_dbsnp_alleles(c("rs1", "rs2")), "must be a single character string")
})

skip_if_no_internet <- function() {
  ok <- tryCatch({
    resp <- httr::GET("https://rest.ensembl.org/", httr::timeout(5))
    httr::status_code(resp) < 500
  }, error = function(e) FALSE)
  if (!ok) testthat::skip("No internet connection / Ensembl REST API unreachable")
}

test_that("fetch_dbsnp_alleles retrieves raw alleles for rs1801133 (MTHFR C677T)", {
  skip_if_no_internet()
  # Confirmed live (2026-07-19): this locus is genuinely triallelic on the
  # Ensembl + strand (G/A/C, ancestral G) -- see test-resolve_pira_alleles.R
  # for the biological explanation. fetch_dbsnp_alleles() is a raw pass-
  # through and must NOT apply any strand conversion itself.
  result <- fetch_dbsnp_alleles("rs1801133")
  expect_setequal(result$alleles, c("G", "A", "C"))
  expect_equal(result$ancestral_allele, "G")
  expect_equal(result$var_class, "SNP")
  expect_equal(result$rsid, "rs1801133")
  expect_equal(result$chromosome, "1")
})

test_that("fetch_dbsnp_alleles errors clearly on an unknown/malformed rsid", {
  skip_if_no_internet()
  expect_error(fetch_dbsnp_alleles("rs999999999999"), "not found")
})
