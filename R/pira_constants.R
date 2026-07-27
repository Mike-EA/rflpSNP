#' PIRA-PCR design constants
#'
#' Default parameters for PIRA-PCR (mismatch primer-introduced restriction
#' analysis) primer design, shared across the PIRA-PCR module functions.
#' Defined here (Phase 2) for reuse starting in Phase 3
#' (\code{scan_pira_candidates()}).
#'
#' \code{PIRA_MISMATCH_WINDOW} gives the allowed range of positions, counted
#' from the primer's 3' end (position 1 = terminal base), where an
#' artificial mismatch may be introduced to create or abolish a restriction
#' site. Default: positions 1 to 3, the range most commonly used in the
#' PIRA-PCR / mismatch-PCR-RFLP literature.
#'
#' @examples
#' PIRA_MISMATCH_WINDOW
#'
#' @export
PIRA_MISMATCH_WINDOW <- 1:3
