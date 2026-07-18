#' Read a gene sequence from a FASTA file
#'
#' Loads a reference (DNA) sequence from a FASTA file and returns it as a
#' `Biostrings::DNAString` object, ready for the rest of the `rflpSNP`
#' functions.
#'
#' @param fasta_file Path to the FASTA file (`.fa`/`.fasta`) containing the
#'   reference sequence (e.g. the full sequence of a gene downloaded from
#'   NCBI).
#' @param seq_index Integer. Index of the sequence to extract if the file
#'   contains more than one record (default `1`).
#'
#' @return A `Biostrings::DNAString` object with the loaded sequence.
#'
#' @details
#' If the FASTA file contains multiple sequences, the function emits an
#' informative message and uses only the sequence indicated by `seq_index`.
#' The input FASTA is expected to contain a single reference sequence
#' (e.g. the genomic region of interest), since the rest of the workflow
#' (`locate_snp()`, `design_primers()`, `simulate_pcr()`) assumes a single
#' reference sequence per analysis.
#'
#' @examples
#' \dontrun{
#' gene_seq <- read_gene_fasta("MTHFR_completeseq.fa")
#' }
#'
#' @export
read_gene_fasta <- function(fasta_file, seq_index = 1) {
  if (!file.exists(fasta_file)) {
    stop(sprintf("FASTA file not found: '%s'", fasta_file))
  }

  full_sequence <- Biostrings::readDNAStringSet(fasta_file)

  if (length(full_sequence) == 0) {
    stop("The FASTA file does not contain any sequence.")
  }
  if (seq_index > length(full_sequence)) {
    stop(sprintf(
      "The FASTA file contains %d sequence(s); seq_index = %d is out of range.",
      length(full_sequence), seq_index
    ))
  }
  if (length(full_sequence) > 1) {
    message(sprintf(
      "The FASTA file contains %d sequences; using sequence %d ('%s').",
      length(full_sequence), seq_index, names(full_sequence)[seq_index]
    ))
  }

  gene_seq <- full_sequence[[seq_index]]
  message(sprintf("Sequence loaded. Length: %d bp", length(gene_seq)))
  gene_seq
}
