# Generic, recursive extractor of the Tm value returned by
# TmCalculator::tm_calculate(). This is needed because, depending on the
# installed version of the package, the result may come back as a GRanges
# object (with the value in mcols()$Tm), as a list with a "Tm" element, or
# nested inside another list. This function walks the returned structure
# and returns the first numeric Tm value it finds, or NULL if none is
# found.
#' @keywords internal
.find_tm_generic <- function(obj, depth = 0) {
  if (depth > 6 || is.null(obj)) return(NULL)

  # 1) S4 objects with metadata columns (GRanges, S4Vectors DataFrame)
  if (isS4(obj)) {
    mc <- tryCatch(S4Vectors::mcols(obj), error = function(e) NULL)
    if (!is.null(mc) && "Tm" %in% colnames(mc)) {
      val <- suppressWarnings(as.numeric(mc$Tm[1]))
      if (!is.na(val)) return(val)
    }
    val <- tryCatch(suppressWarnings(as.numeric(obj$Tm[1])), error = function(e) NA)
    if (!is.na(val)) return(val)
  }

  # 2) A "Tm" attribute attached to the object
  a <- attr(obj, "Tm")
  if (!is.null(a)) {
    val <- suppressWarnings(as.numeric(a[1]))
    if (!is.na(val)) return(val)
  }

  # 3) Lists: look for an element named "Tm"; if it is not numeric, recurse
  if (is.list(obj)) {
    if (!is.null(names(obj)) && "Tm" %in% names(obj)) {
      candidate <- obj[["Tm"]]
      if (is.numeric(candidate) && length(candidate) >= 1) {
        return(as.numeric(candidate[1]))
      }
      val <- .find_tm_generic(candidate, depth + 1)
      if (!is.null(val)) return(val)
    }
    for (elem in obj) {
      val <- .find_tm_generic(elem, depth + 1)
      if (!is.null(val)) return(val)
    }
  }

  return(NULL)
}

#' Calculate the melting temperature (Tm) of an oligonucleotide
#'
#' Calculates the Tm of a primer using the nearest-neighbor method
#' (Nearest-Neighbor) via the `TmCalculator` package, with salt correction.
#'
#' @param seq_dna Oligonucleotide sequence, as a `DNAString` or character
#'   string.
#' @param Na Na+ concentration in mM (default `100`).
#' @param Mg Mg2+ concentration in mM (default `2`).
#' @param dNTPs Total dNTP concentration in mM (default `0.2`).
#' @param oligo_conc_nM Oligo concentration in nM (default `500`).
#' @param nn_table Nearest-neighbor thermodynamic table to use (default
#'   `"DNA_NN_SantaLucia_2004"`).
#' @param salt_corr_method Salt-correction formula to use (default
#'   `"Owczarzy2004"`, appropriate for PCR with Mg2+ and dNTPs present). See
#'   [compare_tm_conditions()] to compare the effect of different
#'   conditions on the same primer.
#'
#' @return The Tm in degrees Celsius (numeric, rounded to 1 decimal), or
#'   `NA` if the calculation could not be performed.
#'
#' @details
#' # Why the Tm calculated here may differ from SnapGene / IDT OligoAnalyzer
#'
#' It is normal and expected for different tools to report different Tm
#' values for the same primer, even when all of them use the
#' nearest-neighbor method. The most common causes, roughly in order of
#' typical impact:
#'
#' 1. **Different salt conditions.** This package defaults to a realistic
#'    PCR (`Na = 100`, `Mg = 2`, `dNTPs = 0.2`), while many online
#'    calculators default to a condition without Mg2+ (e.g. `Na = 50`,
#'    `Mg = 0`). Mg2+ stabilizes the duplex much more per mM than Na+, so
#'    its mere presence can raise the Tm by several degrees. **This is by
#'    far the most likely reason `rflpSNP` reports a higher Tm**: it is not
#'    necessarily an error, but a more realistic salt condition for a PCR
#'    with a standard buffer (which always contains Mg2+). For a fair
#'    comparison, enter exactly the same Na+, Mg2+, dNTP and oligo
#'    concentration values in both tools.
#' 2. **Assumed oligo concentration.** A more concentrated oligo
#'    (higher `oligo_conc_nM`) predicts a slightly higher Tm. IDT
#'    OligoAnalyzer defaults to 0.25 uM (250 nM); the default here is
#'    500 nM (closer to a primer's concentration in an already-prepared
#'    PCR mix).
#' 3. **Salt-correction formula.** `Owczarzy2004` (used here) and other
#'    formulas (`Owczarzy2008`, `SantaLucia1998`, etc.) do not give exactly
#'    the same result, even with identical concentrations.
#' 4. **Nearest-neighbor thermodynamic parameter table.** SantaLucia 2004,
#'    SantaLucia & Hicks 1998, or Sugimoto 1996 differ slightly in their
#'    enthalpy/entropy values per dinucleotide.
#'
#' `TmCalculator::tm_calculate()` has also changed its interface across
#' package versions (for example, the salt-correction argument has been
#' called `salt_method` in some versions and `salt_corr_method` in others,
#' and the output object has moved from simple lists to `GRanges` objects).
#' This function first tries the argument name current on CRAN
#' (`salt_corr_method`) and, if that fails, retries with the name used in
#' earlier versions (`salt_method`) before extracting the Tm value with a
#' generic extractor that tolerates both output formats. Run
#' [check_tm_backend()] once to diagnose your installation if you are
#' unsure which version is active.
#'
#' @seealso [check_tm_backend()], [compare_tm_conditions()]
#' @export
calc_tm <- function(seq_dna, Na = 100, Mg = 2, dNTPs = 0.2, oligo_conc_nM = 500,
                     nn_table = "DNA_NN_SantaLucia_2004",
                     salt_corr_method = "Owczarzy2004") {
  seq_char <- if (inherits(seq_dna, "DNAString")) as.character(seq_dna) else seq_dna
  seq_char <- toupper(seq_char)

  try_call <- function(salt_arg_name) {
    args <- list(
      input_seq = seq_char, method = "tm_nn", nn_table = nn_table,
      dnac_high = oligo_conc_nM, dnac_low = oligo_conc_nM,
      Na = Na, Mg = Mg, dNTPs = dNTPs
    )
    args[[salt_arg_name]] <- salt_corr_method
    do.call(TmCalculator::tm_calculate, args)
  }

  result <- tryCatch(try_call("salt_corr_method"), error = function(e) NULL)
  if (is.null(result)) {
    result <- tryCatch(try_call("salt_method"), error = function(e) {
      message("calc_tm(): tm_calculate() failed with both known salt-correction argument names: ", conditionMessage(e))
      NULL
    })
  }

  if (is.null(result)) return(NA_real_)

  tm_value <- .find_tm_generic(result)
  if (is.null(tm_value)) return(NA_real_)
  round(tm_value, 1)
}

#' Compare the Tm of a primer under different condition profiles
#'
#' Calculates the Tm of the same primer under several typical salt/oligo
#' condition profiles, to help understand how much of a discrepancy with
#' another tool (SnapGene, IDT OligoAnalyzer, etc.) is explained by the
#' assumed conditions rather than by the method itself.
#'
#' @param seq_dna Primer sequence (character or `DNAString`).
#' @param profiles Named list of profiles to compare. Each profile is a list
#'   with `Na`, `Mg`, `dNTPs` and `oligo_conc_nM`. If `NULL` (default),
#'   three reference profiles are used (see Details).
#' @param salt_corr_method Salt-correction formula to use for all profiles
#'   (default `"Owczarzy2004"`).
#'
#' @details
#' Default profiles included:
#' \describe{
#'   \item{`typical_PCR`}{`Na=100, Mg=2, dNTPs=0.2, oligo=500 nM` - the
#'     default used by [calc_tm()] and [design_primers()]; represents a PCR
#'     with a standard buffer (with Mg2+).}
#'   \item{`no_Mg2`}{`Na=50, Mg=0, dNTPs=0, oligo=250 nM` - approximates the
#'     default condition of several online calculators when no PCR buffer
#'     with Mg2+ is specified.}
#'   \item{`low_salt`}{`Na=50, Mg=0, dNTPs=0, oligo=50 nM`.}
#' }
#' These profiles are approximations for classroom exploration. For an
#' exact numeric comparison against a specific tool, set `Na`, `Mg`,
#' `dNTPs` and `oligo_conc_nM` to match exactly the salt/oligo values
#' configured in that tool (not just the calculation method).
#'
#' @return A `data.frame` with columns `profile`, `Na`, `Mg`, `dNTPs`,
#'   `oligo_nM` and `tm`, one row per profile.
#'
#' @examples
#' \dontrun{
#' compare_tm_conditions("TGGTCTCTTCATCCCTCGCCTTGAA")
#' }
#'
#' @seealso [calc_tm()]
#' @export
compare_tm_conditions <- function(seq_dna, profiles = NULL,
                                   salt_corr_method = "Owczarzy2004") {
  if (is.null(profiles)) {
    profiles <- list(
      typical_PCR = list(Na = 100, Mg = 2, dNTPs = 0.2, oligo_conc_nM = 500),
      no_Mg2      = list(Na = 50,  Mg = 0, dNTPs = 0,   oligo_conc_nM = 250),
      low_salt    = list(Na = 50,  Mg = 0, dNTPs = 0,   oligo_conc_nM = 50)
    )
  }

  rows <- lapply(names(profiles), function(name) {
    p <- profiles[[name]]
    tm <- calc_tm(
      seq_dna, Na = p$Na, Mg = p$Mg, dNTPs = p$dNTPs,
      oligo_conc_nM = p$oligo_conc_nM, salt_corr_method = salt_corr_method
    )
    data.frame(
      profile = name, Na = p$Na, Mg = p$Mg, dNTPs = p$dNTPs,
      oligo_nM = p$oligo_conc_nM, tm = tm, stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' Diagnose the installed version of TmCalculator
#'
#' Runs a test call to `TmCalculator::tm_calculate()` and prints the
#' returned structure. Useful for students installing the package for the
#' first time who want to confirm that [calc_tm()] will be able to extract
#' the Tm correctly on their machine, or for debugging if [calc_tm()]
#' unexpectedly returns `NA`.
#'
#' @return Invisible; prints a diagnostic message to the console.
#' @export
check_tm_backend <- function() {
  rflp_version <- tryCatch(
    as.character(utils::packageVersion("rflpSNP")),
    error = function(e) "development"
  )

  tm_version <- tryCatch(
    as.character(utils::packageVersion("TmCalculator")),
    error = function(e) "not installed"
  )

  tm_arguments <- tryCatch(
    names(formals(TmCalculator::tm_calculate)),
    error = function(e) character()
  )

  salt_argument <- intersect(
    c("salt_corr_method", "salt_method"),
    tm_arguments
  )

  cat("\n=== rflpSNP Tm backend check ===\n")
  cat("R version            :", R.version.string, "\n")
  cat("rflpSNP version      :", rflp_version, "\n")
  cat("TmCalculator version :", tm_version, "\n")

  if (length(salt_argument) > 0) {
    cat("Salt argument        :", salt_argument[1], "\n")
  } else {
    cat("Salt argument        : not recognized\n")
  }

  test_tm <- tryCatch(
    calc_tm("ATGCGATGCGATGCATGCA"),
    error = function(e) {
      cat("Backend error        :", conditionMessage(e), "\n")
      NA_real_
    }
  )

  if (is.numeric(test_tm) &&
      length(test_tm) == 1 &&
      is.finite(test_tm)) {
    cat(sprintf(
      "[OK] calc_tm() is working correctly; test Tm = %.1f\u00B0C\n\n",
      test_tm
    ))
  } else {
    cat("\n[DIAGNOSTIC rflpSNP]\n")
    cat("The Tm backend did not return a valid numeric value.\n")
    cat("Validated TmCalculator version: 1.0.8\n")
    cat("Restart R and reinstall rflpSNP with its required dependencies:\n\n")
    cat(
      'remotes::install_github("Mike-EA/rflpSNP", ',
      'dependencies = NA, upgrade = "never", force = TRUE)\n\n',
      sep = ""
    )
  }

  invisible(test_tm)
}
