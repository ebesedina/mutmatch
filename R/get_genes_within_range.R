#' Get Genes Within a Given Range on a Chromosome
#'
#' This function takes a chromosome name and a genomic range specified by left and right borders.
#' It returns the Ensembl IDs of genes that fall within this range.
#'
#' @param chromosome The name of the chromosome.
#' @param leftBorder The leftmost genomic coordinate.
#' @param rightBorder The rightmost genomic coordinate.
#' @return A vector containing Ensembl IDs of genes within the specified range.
#' @export
#' @examples
#' \dontrun{
#' get_genes_within_range("chr1", 1000000, 2000000)
#' }
get_genes_within_range <- function(chromosome, leftBorder, rightBorder) {
  # Create a GenomicRanges object for the input range
  grObject <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(chromosome),
    ranges = IRanges::IRanges(start = leftBorder, end = rightBorder)
  )

  # Filter for only gene type records from the mutmatch_gtf_grch37 dataset
  gtf_genes <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type == "gene"]

  # Match coordinate styles between datasets
  GenomeInfoDb::seqlevelsStyle(gtf_genes) <- "NCBI"
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(grObject) <- "NCBI")

  # Find overlapping genes ignoring the strand
  genes <- GenomicRanges::findOverlaps(query = grObject, subject = gtf_genes, ignore.strand = TRUE)

  # Convert to data table and filter for protein-coding genes
  genes <- data.table::as.data.table(gtf_genes[S4Vectors::subjectHits(genes)])
  genes <- genes[gene_biotype == "protein_coding"]$gene_id

  # Filter for genes present in the topTranscripts dataset
  genes <- genes[genes %in% mutmatch_topTranscripts$Gene]

  return(genes)
}
