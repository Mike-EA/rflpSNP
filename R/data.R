#' Restriction Enzyme Database
#'
#' A curated dataset containing 610 restriction enzymes, including their
#' recognition sites, cut positions, and supplier availability, adapted for
#' PCR-RFLP and PIRA-PCR design.
#'
#' @format A data frame with 610 rows and 13 variables:
#' \describe{
#'   \item{enzyme}{Enzyme name, e.g., "EcoRI"}
#'   \item{recognition_site}{Recognition sequence in IUPAC format}
#'   \item{site_length}{Length of the recognition site in base pairs}
#'   \item{fst5, fst3}{Cut positions relative to the recognition motif}
#'   \item{ovhg, ovhg_seq}{Overhang length and sequence}
#'   \item{cut_type}{Type of cut, e.g., "blunt", "5_overhang"}
#'   \item{palindromic}{Logical indicating if the site is palindromic}
#'   \item{methylation_sensitive}{Logical indicating methylation sensitivity}
#'   \item{n_suppliers}{Number of commercial suppliers}
#'   \item{suppliers}{List of suppliers}
#'   \item{rebase_id}{REBASE database identifier}
#' }
#' @source REBASE / Commercial Catalogs
"enzymes_db"
