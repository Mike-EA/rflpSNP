# ==============================================================================
# BUILD SCRIPT: restriction_enzymes dataset (rflpSNP, PIRA-PCR module, Phase 1)
# ==============================================================================
# Reads data-raw/restriction_enzymes_source.csv (curated from Biopython's
# Restriction_Dictionary.py, REBASE emboss files v404/2024; see
# data-raw/PROVENANCE.md for full methodology and filtering criteria) and
# builds the package-internal `restriction_enzymes` data object.
#
# Run this script from the package root, e.g.:
#   source("data-raw/build_restriction_enzymes.R")
#
# Requires: usethis (for use_data), otherwise falls back to save().
# ==============================================================================

library(utils)

csv_path <- file.path("data-raw", "restriction_enzymes_source.csv")
stopifnot(file.exists(csv_path))

raw <- read.csv(csv_path, stringsAsFactors = FALSE, colClasses = c(
  enzyme                = "character",
  recognition_site      = "character",
  site_length            = "integer",
  fst5                   = "integer",
  fst3                   = "integer",
  ovhg                   = "integer",
  ovhg_seq               = "character",
  cut_type               = "character",
  palindromic            = "logical",
  methylation_sensitive  = "logical",
  n_suppliers            = "integer",
  suppliers              = "character",
  rebase_id              = "integer"
))

# --- Validation --------------------------------------------------------------

iupac_alphabet <- strsplit("ACGTRYSWKMBDHVN", "")[[1]]

validate_restriction_enzymes <- function(df) {
  problems <- character()

  if (anyDuplicated(df$enzyme) > 0) {
    problems <- c(problems, sprintf(
      "Duplicate enzyme names: %s",
      paste(df$enzyme[duplicated(df$enzyme)], collapse = ", ")
    ))
  }

  bad_chars <- vapply(df$recognition_site, function(s) {
    chars <- strsplit(toupper(s), "")[[1]]
    any(!chars %in% iupac_alphabet)
  }, logical(1))
  if (any(bad_chars)) {
    problems <- c(problems, sprintf(
      "Non-IUPAC characters in recognition_site for: %s",
      paste(df$enzyme[bad_chars], collapse = ", ")
    ))
  }

  len_mismatch <- nchar(df$recognition_site) != df$site_length
  if (any(len_mismatch)) {
    problems <- c(problems, sprintf(
      "site_length mismatch with recognition_site for: %s",
      paste(df$enzyme[len_mismatch], collapse = ", ")
    ))
  }

  bad_cut_type <- !df$cut_type %in% c("blunt", "5_overhang", "3_overhang")
  if (any(bad_cut_type)) {
    problems <- c(problems, sprintf(
      "Unrecognized cut_type value for: %s",
      paste(df$enzyme[bad_cut_type], collapse = ", ")
    ))
  }

  # Cross-check cut_type against ovhg sign (see PROVENANCE.md for convention)
  expected_cut_type <- ifelse(df$ovhg == 0, "blunt",
                        ifelse(df$ovhg > 0, "3_overhang", "5_overhang"))
  mismatch_sign <- expected_cut_type != df$cut_type
  if (any(mismatch_sign)) {
    problems <- c(problems, sprintf(
      "cut_type does not match ovhg sign convention for: %s",
      paste(df$enzyme[mismatch_sign], collapse = ", ")
    ))
  }

  too_short <- df$site_length < 4
  if (any(too_short)) {
    problems <- c(problems, sprintf(
      "Recognition site shorter than 4 bp (unexpected after filtering): %s",
      paste(df$enzyme[too_short], collapse = ", ")
    ))
  }

  if (length(problems) > 0) {
    stop("restriction_enzymes validation failed:\n- ",
         paste(problems, collapse = "\n- "))
  }

  invisible(TRUE)
}

validate_restriction_enzymes(raw)
cat(sprintf("Validation passed: %d enzymes.\n", nrow(raw)))

# --- Reference sanity check against known textbook enzymes -------------------
# Confirms overhang sign convention and recognition sites are as expected.
# See data-raw/PROVENANCE.md for how this was originally established.

reference_checks <- list(
  EcoRI   = list(site = "GAATTC", cut_type = "5_overhang", ovhg_seq = "AATT"),
  BamHI   = list(site = "GGATCC", cut_type = "5_overhang", ovhg_seq = "GATC"),
  HindIII = list(site = "AAGCTT", cut_type = "5_overhang", ovhg_seq = "AGCT"),
  PstI    = list(site = "CTGCAG", cut_type = "3_overhang", ovhg_seq = "TGCA"),
  SmaI    = list(site = "CCCGGG", cut_type = "blunt",       ovhg_seq = ""),
  HinfI   = list(site = "GANTC",  cut_type = "5_overhang", ovhg_seq = "ANT")
)

for (enz in names(reference_checks)) {
  row <- raw[raw$enzyme == enz, ]
  if (nrow(row) != 1) stop("Reference enzyme missing from dataset: ", enz)
  exp <- reference_checks[[enz]]
  if (row$recognition_site != exp$site ||
      row$cut_type != exp$cut_type ||
      row$ovhg_seq != exp$ovhg_seq) {
    stop("Reference enzyme mismatch for ", enz, " - check source data/convention.")
  }
}
cat("Reference enzyme sanity check passed (EcoRI, BamHI, HindIII, PstI, SmaI, HinfI).\n")

# --- Finalize object -----------------------------------------------------

restriction_enzymes <- raw
restriction_enzymes$palindromic <- as.logical(restriction_enzymes$palindromic)
restriction_enzymes$methylation_sensitive <- as.logical(restriction_enzymes$methylation_sensitive)
rownames(restriction_enzymes) <- NULL

# --- Save ------------------------------------------------------------------

if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(restriction_enzymes, overwrite = TRUE)
} else {
  if (!dir.exists("data")) dir.create("data")
  save(restriction_enzymes, file = file.path("data", "restriction_enzymes.rda"),
       compress = "bzip2")
}

cat(sprintf(
  "\nSaved restriction_enzymes: %d enzymes, %d columns -> data/restriction_enzymes.rda\n",
  nrow(restriction_enzymes), ncol(restriction_enzymes)
))
