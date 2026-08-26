# Synthetic fixtures for the ARMS-PCR foundation tests.  They deliberately
# contain no external identifiers or data files, so their coordinates are
# stable and can be checked by hand.

arms_fixture_forward <- function() {
  Biostrings::DNAString("AACCGGTTAACCGGTTACGA")
}

arms_fixture_reverse <- function() {
  # The forward reference contains ACGATC at positions 6--11. Its reverse
  # complement is the flank GATCGT, which is found only on the other strand.
  Biostrings::DNAString("TTGCAACGATCCGGAATTC")
}

arms_fixture_repeated_flank <- function() {
  Biostrings::DNAString("AACCGGTAACCGGTC")
}

arms_fixture_unique_amplicon <- function() {
  Biostrings::DNAString("ATGCGTACCTAGGCTAACGTTGCA")
}
