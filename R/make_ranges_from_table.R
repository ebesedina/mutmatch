#' Create Genomic Ranges Object from Genes Data Table
#'
#' This function takes a data table containing gene information, such as genomic start, genomic end, and chromosome,
#' and returns a GenomicRanges object.
#'
#' @param rangesData A data frame containing columns 'chromosome_name', 'genomic_coding_start', 'genomic_coding_end',
#'        'strand', 'length', 'start_position', and 'end_position'.
#' @param exonsCoordinates A logical flag specifying if exon coordinates should be used. Default is TRUE.
#'
#' @return A GenomicRanges object containing the gene information.
#'
#' @examples
#' \dontrun{
#' genesData <- data.frame(
#'   chromosome_name = c("1", "1"),
#'   genomic_coding_start = c(1000, 2000),
#'   genomic_coding_end = c(1500, 2500),
#'   strand = c("+", "-"),
#'   length = c(500, 500),
#'   start_position = c(1000, 2000),
#'   end_position = c(1500, 2500)
#' )
#' grGenesData <- make_ranges_from_table(genesData)
#' }
#'
#' @export
make_ranges_from_table <- function(rangesData, exonsCoordinates = TRUE) {
  # Create GenomicRanges object from the input data frame
  if (!exonsCoordinates) {
    grObject <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(rangesData$chromosome),
      ranges = IRanges::IRanges(
        start = rangesData$gene_start,
        end = rangesData$gene_end
      )
    ) %>% base::unique()
  } else {
    grObject <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(rangesData$chromosome),
      ranges = IRanges::IRanges(
        start = rangesData$exon_start,
        end = rangesData$exon_end
      ),
      strand = rangesData$strand,
      cds_length = rangesData$length,
      gene_start = rangesData$gene_start,
      gene_end = rangesData$gene_end
    )

    # Extend the start and end coordinates by 5 nucleotides
    BiocGenerics::start(grObject) <- BiocGenerics::start(grObject) - 5
    BiocGenerics::end(grObject) <- BiocGenerics::end(grObject) + 5
  }

  # Standardize naming of the chromosomes
  GenomeInfoDb::seqlevelsStyle(grObject) <- "UCSC"

  # Return the GenomicRanges object
  return(grObject)
}
