# IUPAC ambiguity code table: each code maps to the set of concrete bases
# it can represent.
#' @keywords internal
.IUPAC_MAP <- list(
  A = "A", C = "C", G = "G", T = "T",
  R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"),
  B = c("C", "G", "T"), D = c("A", "G", "T"), H = c("A", "C", "T"), V = c("A", "C", "G"),
  N = c("A", "C", "G", "T")
)

# Returns the concrete bases represented by an IUPAC code.
#' @keywords internal
.iupac_expand <- function(code) {
  code <- toupper(code)
  if (!code %in% names(.IUPAC_MAP)) {
    stop(sprintf("Unrecognized IUPAC code: '%s'", code))
  }
  .IUPAC_MAP[[code]]
}

# Returns TRUE if a concrete base satisfies an IUPAC code.
#' @keywords internal
.iupac_matches <- function(base, code) {
  toupper(base) %in% .iupac_expand(code)
}

# Complement of a single concrete base (A/C/G/T only).
#' @keywords internal
.complement_base <- function(base) {
  map <- c(A = "T", T = "A", C = "G", G = "C")
  unname(map[toupper(base)])
}

# IUPAC-aware complement of a single ambiguity code (e.g. R <-> Y).
#' @keywords internal
.IUPAC_COMPLEMENT <- c(
  A = "T", T = "A", C = "G", G = "C",
  R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
  B = "V", V = "B", D = "H", H = "D", N = "N"
)

#' @keywords internal
.iupac_complement_code <- function(code) unname(.IUPAC_COMPLEMENT[toupper(code)])

# Reverse-complements a motif (character vector of IUPAC codes, 5'->3' in,
# 5'->3' out). For a palindromic motif this returns the same string; for a
# non-palindromic one it returns the sequence as it would read on the
# opposite strand -- required to correctly search the "downstream" (mirrored)
# side in find_pira_sites().
#' @keywords internal
.iupac_revcomp_motif <- function(motif_vec) {
  rev(vapply(motif_vec, .iupac_complement_code, character(1), USE.NAMES = FALSE))
}

#' Load a custom restriction enzyme panel from a CSV file
#'
#' Reads a CSV file of restriction enzymes (e.g. a REBASE-derived curated
#' database) and adapts it to the `name`/`motif`/`cut_offset` schema used
#' by [find_pira_sites()] and [design_pira_primers()].
#'
#' @param path Path to the CSV file.
#' @param name_col,motif_col,cut_offset_col Column names in the CSV holding
#'   the enzyme name, the recognition site (IUPAC), and the cut position
#'   (counted from the first base of the site, 5' end, same convention as
#'   REBASE/Biopython's `fst5`; `0` means the enzyme cuts before the first
#'   base of the site). Defaults match a REBASE/Biopython-style export:
#'   `"enzyme"`, `"recognition_site"`, `"fst5"`.
#' @param commercial_only If `TRUE` (default) and a `n_suppliers` column is
#'   present, keep only enzymes with at least one commercial supplier
#'   (`n_suppliers >= 1`).
#' @param max_site_length Enzymes with a recognition site longer than this
#'   are dropped (default `12`); very long/degenerate sites rarely help
#'   PIRA-PCR design and slow down the search substantially.
#'
#' @return A `data.frame` with columns `name`, `motif`, `cut_offset`, ready
#'   to pass as `enzyme_panel` to [find_pira_sites()] or
#'   [design_pira_primers()].
#'
#' @examples
#' \dontrun{
#' panel <- read_enzyme_panel("enzymes_610.csv")
#' find_pira_sites(gene_seq, snp_pos, "C", "T", enzyme_panel = panel)
#' }
#'
#' @export
read_enzyme_panel <- function(path, name_col = "enzyme", motif_col = "recognition_site",
                               cut_offset_col = "fst5", commercial_only = TRUE,
                               max_site_length = 12) {
  raw <- utils::read.csv(path, stringsAsFactors = FALSE)

  missing_cols <- setdiff(c(name_col, motif_col, cut_offset_col), colnames(raw))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in '%s': %s. Available columns: %s",
      path, paste(missing_cols, collapse = ", "), paste(colnames(raw), collapse = ", ")
    ))
  }

  panel <- data.frame(
    name = raw[[name_col]],
    motif = toupper(raw[[motif_col]]),
    cut_offset = suppressWarnings(as.integer(raw[[cut_offset_col]])),
    stringsAsFactors = FALSE
  )

  n0 <- nrow(panel)
  panel <- panel[!is.na(panel$cut_offset) & nchar(panel$motif) > 0, ]
  panel <- panel[nchar(panel$motif) <= max_site_length, ]

  if (commercial_only && "n_suppliers" %in% colnames(raw)) {
    keep <- raw$n_suppliers[match(panel$name, raw[[name_col]])] >= 1
    panel <- panel[!is.na(keep) & keep, ]
  }

  panel <- panel[!duplicated(panel$name), ]
  rownames(panel) <- NULL

  message(sprintf(
    "Loaded %d enzyme(s) from '%s' (%d dropped: missing/invalid data, site too long, or no commercial supplier).",
    nrow(panel), path, n0 - nrow(panel)
  ))

  panel
}

#' Restriction enzyme panel used by default in PIRA-PCR design
#'
#' Returns the restriction enzyme panel used as the default `enzyme_panel`
#' in [find_pira_sites()] and [design_pira_primers()].
#'
#' @param source `"full"` (default) loads the curated panel bundled with
#'   the package (`inst/extdata/restriction_enzymes.csv`, several hundred
#'   enzymes). `"builtin"` returns a small, hand-picked panel of ~33 common
#'   commercially available enzymes, useful for quick examples/teaching or
#'   as a fallback if the bundled file is unavailable for some reason
#'   (e.g. a development install missing `inst/extdata/`).
#' @param commercial_only,max_site_length Passed to [read_enzyme_panel()]
#'   when `source = "full"` (see that function for defaults/meaning).
#'
#' @details
#' The bundled panel (`source = "full"`) is expected to be REBASE-derived;
#' see `inst/NOTICE.md` for attribution and licensing (REBASE data is
#' distributed under CC BY-NC -- non-commercial use). You can always
#' extend or replace either panel: build your own `data.frame` with the
#' same three columns (`name`, `motif`, `cut_offset`) -- `motif` uses IUPAC
#' codes (`N`, `R`, `Y`, etc.); `cut_offset` is the position of the cut
#' within the motif, counted from its first base (1-based), matching the
#' convention used by [find_restriction_site()] -- and pass it directly as
#' `enzyme_panel` to [find_pira_sites()]/[design_pira_primers()], bypassing
#' this function entirely.
#'
#' @return A `data.frame` with columns `name`, `motif`, `cut_offset`.
#' @export
pira_enzyme_panel <- function(source = c("full", "builtin"),
                               commercial_only = TRUE, max_site_length = 12) {
  source <- match.arg(source)

  if (source == "builtin") {
    return(.pira_enzyme_panel_builtin())
  }

  path <- system.file("extdata", "restriction_enzymes.csv", package = "rflpSNP")
  if (!nzchar(path)) {
    warning(
      "Bundled enzyme database (inst/extdata/restriction_enzymes.csv) not found ",
      "in this installation of rflpSNP; falling back to the small builtin panel. ",
      "Use pira_enzyme_panel('builtin') to silence this warning, or reinstall the ",
      "package with the data file in place."
    )
    return(.pira_enzyme_panel_builtin())
  }

  read_enzyme_panel(path, commercial_only = commercial_only, max_site_length = max_site_length)
}

# The original small, hand-curated panel (~33 common enzymes), kept as an
# internal fallback and as the explicit "builtin" option.
#' @keywords internal
.pira_enzyme_panel_builtin <- function() {
  data.frame(
    name = c(
      "EcoRI", "BamHI", "HindIII", "PstI", "SalI", "XhoI", "NotI", "SmaI",
      "KpnI", "SacI", "NcoI", "NdeI", "SphI", "XbaI", "ApaI", "HinfI",
      "MspI", "TaqI", "AluI", "RsaI", "HaeIII", "DdeI", "MboI", "NlaIII",
      "EcoRV", "ScaI", "PvuII", "BglII", "ClaI", "NheI", "SpeI", "StuI", "MseI"
    ),
    motif = c(
      "GAATTC", "GGATCC", "AAGCTT", "CTGCAG", "GTCGAC", "CTCGAG", "GCGGCCGC", "CCCGGG",
      "GGTACC", "GAGCTC", "CCATGG", "CATATG", "GCATGC", "TCTAGA", "GGGCCC", "GANTC",
      "CCGG", "TCGA", "AGCT", "GTAC", "GGCC", "CTNAG", "GATC", "CATG",
      "GATATC", "AGTACT", "CAGCTG", "AGATCT", "ATCGAT", "GCTAGC", "ACTAGT", "AGGCCT", "TTAA"
    ),
    cut_offset = c(
      1, 1, 1, 5, 1, 1, 2, 3,
      5, 5, 1, 2, 5, 1, 5, 1,
      1, 1, 2, 2, 2, 1, 0, 4,
      3, 3, 3, 1, 2, 1, 1, 3, 1
    ),
    stringsAsFactors = FALSE
  )
}
