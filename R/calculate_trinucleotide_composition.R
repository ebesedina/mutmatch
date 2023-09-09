#' Calculate Trinucleotide Composition in a Genomic Region
#'
#' This function computes the trinucleotide composition for a given genomic range object.
#' Both sense and antisense counts are considered.
#'
#' @param grObject A GenomicRanges object representing the genomic region for which trinucleotide frequency will be calculated.
#'
#' @return A data frame containing the trinucleotide contexts and their corresponding counts.
#' @export
#'
#' @examples
#' \dontrun{
#' gr <- GRanges("chr1", IRanges(1000, 2000))
#' result <- calculate_trinucleotide_composition(gr)
#' }
calculate_trinucleotide_composition <- function(grObject) {
  # Adjusting start and end positions of the genomic range
  GenomicRanges::start(grObject) <- GenomicRanges::start(grObject) - 1
  GenomicRanges::end(grObject) <- GenomicRanges::end(grObject) + 1

  # Get the sequence and its reverse complement
  seq_sense <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, grObject) %>% suppressWarnings()
  seq_antisense <- Biostrings::reverseComplement(seq_sense) %>% suppressWarnings()

  # Calculate the trinucleotide frequency
  sense_counts <- Biostrings::trinucleotideFrequency(seq_sense) %>%
    base::as.data.frame() %>%
    base::colSums()
  antisense_counts <- Biostrings::trinucleotideFrequency(seq_antisense) %>%
    base::as.data.frame() %>%
    base::colSums()

  # Consolidate the data into a data frame
  Composition3nt <- data.table::data.table(
    Context = base::names(sense_counts),
    Sense_counts = base::as.numeric(sense_counts),
    Antisense_counts = base::as.numeric(antisense_counts)
  ) %>%
    dplyr::mutate(
      Counts = Sense_counts + Antisense_counts
    ) %>%
    dplyr::filter(
      !(base::substr(Context, 2, 2) %in% base::c("A", "G")) &
        Counts != 0
    )

  return(Composition3nt)
}
