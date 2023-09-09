#' Identify Closest Protein-Coding Neighboring Gene
#'
#' This function searches for the nearest protein-coding gene adjacent to a specified target gene.
#' The search can be conducted either upstream ('left') or downstream ('right') relative to the target,
#' without consideration for gene directionality or orientation on the DNA strand.
#'
#' @param grGene A GenomicRanges object representing the genomic coordinates of the target gene for which to find the closest neighbor.
#' @param grGenes A GenomicRanges object representing the genomic coordinates of the available set of genes within which to conduct the search.
#' @param side A character string specifying the direction of the search: 'left' for upstream or 'right' for downstream. The search is performed without taking into account the orientation of the gene on the DNA strand.
#'
#' @return Returns a vector of unique gene IDs representing the closest protein-coding neighbor(s) to the target gene. Returns NULL if no such genes are found.
#'
#' @examples
#' # Assuming `grTarget` and `grGeneSet` are defined GenomicRanges objects
#' \dontrun{
#' find_closest_neighbor(grTarget, grGeneSet, "left")
#' find_closest_neighbor(grTarget, grGeneSet, "right")
#' }
#' @export
find_closest_neighbor <- function(grGene,
                                  grGenes,
                                  side = c("left", "right")) {
  # Ensure 'side' parameter contains a valid value
  side <- match.arg(side)

  # Match styles
  GenomeInfoDb::seqlevelsStyle(grGene) <- "NCBI"
  GenomeInfoDb::seqlevelsStyle(grGene) <- "NCBI"

  # Search for the closest gene on the left side of the target gene
  if (side == "left") {
    left_genes <- GenomicRanges::follow(
      x = grGene,
      subject = grGenes,
      ignore.strand = TRUE
    )
    if (is.na(left_genes)) {
      return(NULL)
    } else {
      # Convert to data table and filter for 'protein_coding' genes
      left_genes <- data.table::as.data.table(grGenes[left_genes])
      left_genes <- left_genes[gene_biotype == "protein_coding"]$gene_id %>% base::unique()

      return(left_genes)
    }

    # Search for the closest gene on the right side of the target gene
  } else {
    right_genes <- GenomicRanges::precede(
      x = grGene,
      subject = grGenes,
      ignore.strand = TRUE
    )
    if (is.na(right_genes)) {
      return(NULL)
    } else {
      # Convert to data table and filter for 'protein_coding' genes
      right_genes <- data.table::as.data.table(grGenes[right_genes])
      right_genes <- right_genes[gene_biotype == "protein_coding"]$gene_id %>% base::unique()

      return(right_genes)
    }
  }
}
