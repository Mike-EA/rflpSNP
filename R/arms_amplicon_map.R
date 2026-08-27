#' Export nucleotide maps of simulated ARMS-PCR products
#'
#' Writes a reproducible, base-by-base map for every expected ARMS-PCR product.
#' Each map contains the forward strand (5' to 3'), the aligned complementary
#' strand (3' to 5'), FASTA coordinates, SNP position and the primers involved.
#'
#' @param arms_pcr A `rflp_arms_pcr` object from [simulate_arms_pcr()].
#' @param output_file Path of the text map to create.
#' @param line_width Number of nucleotides per displayed row.
#'
#' @return Invisibly, the path of the written map.
#' @export
export_arms_amplicon_map_txt <- function(arms_pcr,
                                         output_file = "arms_pcr_amplicon_map.txt",
                                         line_width = 60L) {
  if (!inherits(arms_pcr, "rflp_arms_pcr")) stop("'arms_pcr' must be returned by simulate_arms_pcr().")
  if (length(line_width) != 1L || is.na(line_width) || line_width < 1L) stop("'line_width' must be a positive integer.")
  set <- arms_pcr$primer_set
  primer_rows <- data.frame(
    name = c("outer_forward", "outer_reverse", "ref_inner", "alt_inner"),
    sequence = c(set$outer_forward_seq, set$outer_reverse_seq, set$ref_inner_seq, set$alt_inner_seq),
    strand = c("forward", "reverse", set$ref_inner_strand, set$alt_inner_strand),
    start = c(set$outer_forward_start, set$outer_reverse_start, set$ref_inner_start, set$alt_inner_start),
    end = c(set$outer_forward_end, set$outer_reverse_end, set$ref_inner_end, set$alt_inner_end),
    mismatch = c(NA_character_, NA_character_, set$ref_mismatch_3prime, set$alt_mismatch_3prime),
    stringsAsFactors = FALSE
  )
  product_primers <- list(
    control = c("outer_forward", "outer_reverse"),
    reference = c("outer_reverse", "ref_inner"),
    alternative = c("outer_forward", "alt_inner")
  )
  report <- c(
    "===========================================================",
    "           TETRA-PRIMER ARMS-PCR NUCLEOTIDE MAP",
    "===========================================================",
    sprintf("SNP coordinate: %d | ref: %s | alt: %s", arms_pcr$snp_pos, arms_pcr$ref_allele, arms_pcr$alt_allele),
    "All genomic coordinates are 1-based and inclusive.",
    "The heterozygote control band represents both allele templates; the displayed sequence is the reference-template representation.",
    ""
  )
  for (genotype in names(arms_pcr$products)) {
    report <- c(report, paste0("--- GENOTYPE: ", genotype, " ---"))
    products <- arms_pcr$products[[genotype]]
    for (i in seq_len(nrow(products))) {
      product <- products[i, ]
      selected <- primer_rows[primer_rows$name %in% product_primers[[product$product]], , drop = FALSE]
      report <- c(report,
        sprintf("Product: %s | template: %s | FASTA: %d-%d | size: %d bp", product$product,
                product$template_allele, product$start, product$end, product$size_bp),
        "Primers (stored 5'->3'):",
        utils::capture.output(print(selected, row.names = FALSE, right = FALSE))
      )
      sequence <- product$sequence
      complement <- as.character(Biostrings::complement(Biostrings::DNAString(sequence)))
      starts <- seq.int(1L, nchar(sequence), by = as.integer(line_width))
      for (offset in starts) {
        last <- min(nchar(sequence), offset + as.integer(line_width) - 1L)
        coordinates <- seq.int(product$start + offset - 1L, product$start + last - 1L)
        snp_marker <- if (arms_pcr$snp_pos %in% coordinates) {
          paste0(strrep(" ", arms_pcr$snp_pos - coordinates[1]), "^")
        } else ""
        report <- c(report,
          sprintf("Position  %d%s%d", coordinates[1], if (length(coordinates) > 1L) "..." else "", coordinates[length(coordinates)]),
          paste0("Forward 5'  ", substr(sequence, offset, last), " 3'"),
          paste0("Complement3'  ", substr(complement, offset, last), " 5'"),
          if (nzchar(snp_marker)) paste0("SNP         ", snp_marker) else ""
        )
      }
      report <- c(report, "")
    }
  }
  writeLines(report, con = output_file)
  message(sprintf("ARMS-PCR nucleotide map exported to: %s", output_file))
  invisible(output_file)
}
