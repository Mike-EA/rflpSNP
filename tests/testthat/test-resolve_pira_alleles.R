test_that("resolve_pira_alleles validates snp_base without network access", {
  expect_error(resolve_pira_alleles(snp_base = "X", rsid = "rs1801133"),
               "must be one of A, C, G, T")
  expect_error(resolve_pira_alleles(snp_base = "CG", rsid = "rs1801133"),
               "must be one of A, C, G, T")
  expect_error(resolve_pira_alleles(snp_base = NULL, rsid = "rs1801133"))
})

skip_if_no_internet <- function() {
  ok <- tryCatch({
    resp <- httr::GET("https://rest.ensembl.org/", httr::timeout(5))
    httr::status_code(resp) < 500
  }, error = function(e) FALSE)
  if (!ok) testthat::skip("No internet connection / Ensembl REST API unreachable")
}

test_that("resolve_pira_alleles correctly resolves MTHFR C677T via the ancestral allele", {
  skip_if_no_internet()
  # Confirmed live (2026-07-19): Ensembl reports rs1801133 as triallelic on
  # the + strand: allele_string "G/A/C", ancestral_allele "G", strand 1.
  # This reflects two distinct real clinical substitutions at the same
  # nucleotide: c.665C>T (classic C677T / A222V) and c.665C>G (A222G).
  # snp_base "C" (the coding-strand/literature reference) does NOT equal the
  # ancestral allele "G" directly, but its complement does -> orientation
  # must resolve to "complement", confirmed via the ancestral allele (not by
  # coincidental direct match, which would silently give the wrong answer
  # here since "C" is also a literal, but unrelated, allele in the raw set).
  res <- resolve_pira_alleles(snp_base = "C", rsid = "rs1801133")

  expect_equal(res$reference_allele, "C")
  expect_equal(res$strand_used, "complement")
  expect_true(res$orientation_verified)
  expect_true(res$multiallelic)
  expect_true("T" %in% res$alternate_alleles)  # c.665C>T, classic variant
  expect_true("G" %in% res$alternate_alleles)  # c.665C>G, rarer variant
  expect_length(res$alternate_alleles, 2)
})
