#' Fetch SNP allele data from Ensembl by rsID
#'
#' Queries the Ensembl REST API (\url{https://rest.ensembl.org}) for a given
#' RefSNP identifier and returns the observed alleles, chromosome position,
#' and related metadata. Used by \code{\link{resolve_pira_alleles}} to obtain
#' the alternate allele needed for PIRA-PCR primer design.
#'
#' This function intentionally avoids depending on the \code{rsnps} package,
#' which has been removed from CRAN. It calls the Ensembl REST API directly
#' via \code{httr} and \code{jsonlite}, both standard, actively maintained
#' CRAN packages.
#'
#' @param rsid Character. A RefSNP identifier, e.g. \code{"rs1801133"}.
#'   Must start with \code{"rs"}.
#' @param species Character. Species alias understood by the Ensembl REST
#'   API. Default \code{"human"}.
#' @param timeout_sec Numeric. Request timeout in seconds. Default 10.
#'
#' @return A list with elements:
#' \describe{
#'   \item{rsid}{Character. The queried/resolved rsID.}
#'   \item{chromosome}{Character. Chromosome name (Ensembl \code{seq_region_name}).}
#'   \item{position}{Integer. Genomic start position of the variant (Ensembl coordinates).}
#'   \item{alleles}{Character vector. All observed alleles at this site (may be > 2).}
#'   \item{ancestral_allele}{Character or \code{NA}. Ancestral allele, if known.}
#'   \item{minor_allele}{Character or \code{NA}. Reported minor allele.}
#'   \item{var_class}{Character. Ensembl variant class, e.g. \code{"SNP"}.}
#'   \item{assembly}{Character. Genome assembly name (e.g. \code{"GRCh38"}).}
#' }
#'
#' @details
#' Errors are reported with informative messages rather than raw \code{httr}
#' errors, covering: no internet connection, rsID not found (HTTP 404), and
#' rate limiting (HTTP 429). If the Ensembl response includes more than one
#' genomic mapping (e.g. alternate scaffolds/patches), the mapping with
#' \code{coord_system == "chromosome"} is used; if none is found, the first
#' mapping is used with a warning.
#'
#' @examples
#' \dontrun{
#' fetch_dbsnp_alleles("rs1801133")  # MTHFR C677T
#' }
#'
#' @export
fetch_dbsnp_alleles <- function(rsid, species = "human", timeout_sec = 10) {

  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required for fetch_dbsnp_alleles(). ",
         "Install it with install.packages(\"httr\").", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for fetch_dbsnp_alleles(). ",
         "Install it with install.packages(\"jsonlite\").", call. = FALSE)
  }

  if (missing(rsid) || is.null(rsid) || !is.character(rsid) || length(rsid) != 1) {
    stop("'rsid' must be a single character string, e.g. \"rs1801133\".", call. = FALSE)
  }
  if (!grepl("^rs[0-9]+$", rsid)) {
    stop("'rsid' does not look like a valid RefSNP identifier: '", rsid,
         "'. Expected a format like \"rs1801133\".", call. = FALSE)
  }

  url <- sprintf(
    "https://rest.ensembl.org/variation/%s/%s?content-type=application/json",
    utils::URLencode(species), utils::URLencode(rsid)
  )

  resp <- tryCatch(
    httr::GET(url, httr::accept("application/json"),
              httr::timeout(timeout_sec)),
    error = function(e) {
      stop("Could not reach the Ensembl REST API. This usually means there ",
           "is no internet connection available. Original error: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  status <- httr::status_code(resp)

  if (status %in% c(400, 404)) {
    stop("rsID '", rsid, "' was not found in Ensembl/dbSNP (HTTP ", status,
         "). Verify the ID, or note that some rsIDs have been merged into ",
         "a different identifier over time.", call. = FALSE)
  }
  if (status == 429) {
    stop("Ensembl REST API rate limit reached (HTTP 429). Wait a moment ",
         "and try again.", call. = FALSE)
  }
  if (status != 200) {
    stop("Ensembl REST API returned an unexpected HTTP status (", status,
         ") for rsID '", rsid, "'.", call. = FALSE)
  }

  body_text <- httr::content(resp, as = "text", encoding = "UTF-8")
  parsed <- tryCatch(
    jsonlite::fromJSON(body_text, simplifyVector = FALSE),
    error = function(e) {
      stop("Failed to parse the Ensembl REST API response as JSON: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  if (is.null(parsed$mappings) || length(parsed$mappings) == 0) {
    stop("Ensembl returned no genomic mappings for rsID '", rsid,
         "'. The variant may not be mapped to the current assembly.",
         call. = FALSE)
  }

  mappings <- parsed$mappings
  coord_systems <- vapply(mappings, function(m) {
    cs <- m$coord_system
    if (is.null(cs)) NA_character_ else cs
  }, character(1))

  chrom_idx <- which(coord_systems == "chromosome")
  if (length(chrom_idx) == 0) {
    warning("No mapping with coord_system == 'chromosome' found for '", rsid,
             "'; using the first available mapping instead. Verify the ",
             "result carefully.", call. = FALSE)
    chosen <- mappings[[1]]
  } else {
    chosen <- mappings[[chrom_idx[1]]]
  }

  allele_string <- chosen$allele_string
  if (is.null(allele_string) || !nzchar(allele_string)) {
    stop("Ensembl mapping for '", rsid, "' has no allele_string field.",
         call. = FALSE)
  }
  alleles <- strsplit(allele_string, "/", fixed = TRUE)[[1]]
  alleles <- toupper(trimws(alleles))

  var_class <- parsed$var_class
  if (is.null(var_class)) var_class <- NA_character_

  if (!identical(var_class, "SNP")) {
    warning("rsID '", rsid, "' has Ensembl var_class = '", var_class,
             "', not 'SNP'. PIRA-PCR primer design assumes a single-",
             "nucleotide variant; results for indels or other variant ",
             "classes are not validated.", call. = FALSE)
  }

  list(
    rsid              = if (is.null(parsed$name)) rsid else parsed$name,
    chromosome        = chosen$seq_region_name %||% NA_character_,
    position          = as.integer(chosen$start %||% NA_integer_),
    alleles           = alleles,
    ancestral_allele  = toupper(chosen$ancestral_allele %||% NA_character_),
    minor_allele      = toupper(parsed$minor_allele %||% NA_character_),
    var_class         = var_class,
    assembly          = chosen$assembly_name %||% NA_character_
  )
}
