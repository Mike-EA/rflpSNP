#' @keywords internal
.arms_extract_set <- function(primer_set) {
  if (inherits(primer_set, "rflp_arms_primers")) primer_set <- primer_set$best_set
  if (!is.data.frame(primer_set) || nrow(primer_set) != 1L) {
    stop("'primer_set' must be one row from design_arms_primers() or a design object with a best set.")
  }
  required <- c("outer_forward_seq", "outer_forward_start", "outer_forward_end",
    "outer_reverse_seq", "outer_reverse_start", "outer_reverse_end", "ref_inner_seq",
    "ref_inner_start", "ref_inner_end", "ref_inner_strand", "ref_mismatch_3prime",
    "alt_inner_seq", "alt_inner_start", "alt_inner_end", "alt_inner_strand", "alt_mismatch_3prime")
  if (!all(required %in% names(primer_set))) stop("'primer_set' does not contain the required ARMS-PCR columns.")
  primer_set
}

#' @keywords internal
.arms_mismatch_index <- function(label, sequence) {
  distance <- suppressWarnings(as.integer(sub("^3'-", "", label)))
  index <- nchar(sequence) - distance
  if (!distance %in% c(2L, 3L) || index < 1L) stop("Invalid deliberate mismatch annotation in 'primer_set'.")
  index
}

#' @keywords internal
.arms_validate_simulation_set <- function(gene_seq, set, snp_pos, ref_allele, alt_allele) {
  input <- .arms_validate_inputs(gene_seq, snp_pos, ref_allele, alt_allele)
  verify_outer <- function(sequence, start, end, strand) {
    if (start < 1L || end > length(gene_seq) || start > end) return(FALSE)
    template <- as.character(Biostrings::subseq(gene_seq, start, end))
    expected <- if (identical(strand, "forward")) template else
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(template)))
    identical(sequence, expected)
  }
  if (!verify_outer(set$outer_forward_seq, set$outer_forward_start, set$outer_forward_end, "forward") ||
      !verify_outer(set$outer_reverse_seq, set$outer_reverse_start, set$outer_reverse_end, "reverse") ||
      set$outer_forward_end >= input$snp_pos || set$outer_reverse_start <= input$snp_pos) {
    stop("Outer-primer sequences, orientation or coordinates are inconsistent with 'gene_seq'.")
  }
  inner <- function(prefix, target, other) {
    mismatch_column <- paste0(sub("_inner$", "", prefix), "_mismatch_3prime")
    primer <- data.frame(sequence = set[[paste0(prefix, "_seq")]],
      genomic_start = set[[paste0(prefix, "_start")]], genomic_end = set[[paste0(prefix, "_end")]],
      strand = set[[paste0(prefix, "_strand")]], mismatch_index = .arms_mismatch_index(
        set[[mismatch_column]], set[[paste0(prefix, "_seq")]]), stringsAsFactors = FALSE)
    endpoint_ok <- if (identical(primer$strand, "forward")) primer$genomic_end == input$snp_pos else
      identical(primer$strand, "reverse") && primer$genomic_start == input$snp_pos
    endpoint_ok && .arms_inner_is_specific(primer, gene_seq, input$snp_pos, target, other)
  }
  if (!inner("ref_inner", input$ref_allele, input$alt_allele) ||
      !inner("alt_inner", input$alt_allele, input$ref_allele)) {
    stop("Inner-primer orientation, coordinates or 3' allele specificity are inconsistent with 'gene_seq'.")
  }
  input
}

#' Simulate tetra-primer ARMS-PCR products for all three genotypes
#'
#' Validates the full four-primer geometry and the intentional terminal
#' mismatch pattern before deriving the expected control and diagnostic
#' products. It does not use the generic mismatch tolerance of [simulate_pcr()].
#'
#' @param gene_seq Reference `Biostrings::DNAString`.
#' @param primer_set One-row set from [design_arms_primers()] or its design object.
#' @param snp_pos One-based SNP coordinate.
#' @param ref_allele,alt_allele Explicit reference and alternative bases.
#'
#' @return A `rflp_arms_pcr` object containing product tables for `ref/ref`,
#'   `ref/alt` and `alt/alt`.
#' @export
simulate_arms_pcr <- function(gene_seq, primer_set, snp_pos, ref_allele, alt_allele) {
  set <- .arms_extract_set(primer_set)
  input <- .arms_validate_simulation_set(gene_seq, set, snp_pos, ref_allele, alt_allele)
  alt_sequence <- Biostrings::DNAString(.arms_replace_base(as.character(gene_seq), input$snp_pos, input$alt_allele))
  product <- function(name, start, end, template, template_allele) data.frame(
    product = name, start = start, end = end, size_bp = end - start + 1L,
    template_allele = template_allele, sequence = as.character(Biostrings::subseq(template, start, end)), stringsAsFactors = FALSE)
  control_ref <- product("control", set$outer_forward_start, set$outer_reverse_end, gene_seq, "ref")
  control_alt <- product("control", set$outer_forward_start, set$outer_reverse_end, alt_sequence, "alt")
  control_het <- control_ref; control_het$template_allele <- "ref/alt"
  allele_product <- function(prefix, name, template, template_allele) {
    if (identical(set[[paste0(prefix, "_strand")]], "forward"))
      product(name, set[[paste0(prefix, "_start")]], set$outer_reverse_end, template, template_allele)
    else product(name, set$outer_forward_start, set[[paste0(prefix, "_end")]], template, template_allele)
  }
  reference <- allele_product("ref_inner", "reference", gene_seq, "ref")
  alternative <- allele_product("alt_inner", "alternative", alt_sequence, "alt")
  result <- list(
    products = list("ref/ref" = rbind(control_ref, reference), "ref/alt" = rbind(control_het, reference, alternative),
                    "alt/alt" = rbind(control_alt, alternative)),
    primer_set = set, snp_pos = input$snp_pos, ref_allele = input$ref_allele, alt_allele = input$alt_allele
  )
  class(result) <- "rflp_arms_pcr"
  result
}

#' Simulate an agarose gel for tetra-primer ARMS-PCR
#'
#' @param arms_pcr Result from [simulate_arms_pcr()].
#' @param ladder_sizes Marker sizes in base pairs.
#'
#' @return A `ggplot` object with marker and the three expected genotype lanes.
#' @export
simulate_arms_gel <- function(arms_pcr, ladder_sizes = c(50, 100, 150, 200, 250, 300, 400, 500)) {
  if (!inherits(arms_pcr, "rflp_arms_pcr")) stop("'arms_pcr' must be returned by simulate_arms_pcr().")
  if (!is.numeric(ladder_sizes) || !length(ladder_sizes) || any(ladder_sizes <= 0)) stop("'ladder_sizes' must contain positive sizes.")
  genotype_order <- c("ref/ref", "ref/alt", "alt/alt")
  samples <- do.call(rbind, lapply(seq_along(genotype_order), function(i) {
    bands <- arms_pcr$products[[genotype_order[i]]]
    data.frame(lane = i + 1L, genotype = genotype_order[i], product = bands$product, size_bp = bands$size_bp)
  }))
  marker <- data.frame(lane = 1L, genotype = "Marker", product = "marker", size_bp = ladder_sizes)
  bands <- rbind(marker, samples); bands$y <- log10(bands$size_bp)
  max_size <- max(bands$size_bp); min_size <- min(bands$size_bp)
  ggplot2::ggplot(bands) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0.5, xmax = 4.5, ymin = log10(min_size * 0.7), ymax = log10(max_size * 1.35)),
      fill = "#0d0f12", color = "grey30") +
    ggplot2::geom_segment(ggplot2::aes(x = lane - .25, xend = lane + .25, y = y, yend = y, colour = product), linewidth = 3.5, lineend = "round") +
    ggplot2::geom_text(data = subset(bands, product != "marker"), ggplot2::aes(x = lane, y = y + .025, label = paste0(size_bp, " bp")), color = "white", size = 3) +
    ggplot2::scale_colour_manual(values = c(marker = "#d8d8d8", control = "#90ee90", reference = "#78b7ff", alternative = "#ffb36b")) +
    ggplot2::scale_x_continuous(breaks = 1:4, labels = c("MW\nmarker", "ref/ref", "ref/alt", "alt/alt")) +
    ggplot2::labs(title = "Virtual tetra-primer ARMS-PCR gel", subtitle = "Control and allele-specific products", x = "Loading lane", y = "Electrophoretic migration (- -> +)", colour = "Band") +
    ggplot2::theme_minimal() + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"), axis.text.x = ggplot2::element_text(face = "bold"))
}
