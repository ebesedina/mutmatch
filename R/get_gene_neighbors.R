#' Identify Neighboring Genes
#'
#' This function identifies genes neighboring a specified target gene, confined to a defined genomic window.
#' It filters these neighboring genes based on a similarity threshold for mutation rates relative to the
#' target gene. In cases where no neighbors are found within the specified window, the function selects
#' the closest upstream and downstream genes that meet the mutation rate similarity criteria.
#'
#' @param grGene The genomic coordinates of the target gene.
#' @param hgnc The HGNC symbol of the target gene.
#' @param neighborsWindow The genomic window within which to search for neighboring genes. Accepts human-readable format like "0.5 Mb".
#' @param outlierNeighborsThreshold The maximum absolute value of the log odds ratio of baseline mutation rate between
#' neighboring genes and the central gene. Genes exceeding this threshold will be filtered out when analyzing neighboring genes.
#'
#' @return A data.table of genomic coordinates for the neighboring genes, or NULL if none are found.
#' @export
get_gene_neighbors <- function(grGene,
                               hgnc,
                               neighborsWindow = "0.5Mb",
                               outlierNeighborsThreshold = 0.2) {
  # Convert the neighborsWindow parameter to actual numerical values (in bp)
  neighborsWindow <- convert_dna_length(neighborsWindow)

  # Extract the chromosome name and the start and end coordinates of the target gene
  chromosome <- GenomeInfoDb::seqnames(grGene) %>% base::as.character()
  start_gene <- BiocGenerics::start(grGene)
  end_gene <- BiocGenerics::end(grGene)

  # Define the search window for neighboring genes
  leftmost_coordinate <- start_gene - neighborsWindow
  left_targetGene_border <- start_gene - 1
  rightmost_coordinate <- end_gene + neighborsWindow
  right_targetGene_border <- end_gene + 1

  # Get genes with similar mutation rates as the target gene
  genesSimMutRate <- get_similar_genes_rmd(hgnc = hgnc, outlierNeighborsThreshold = outlierNeighborsThreshold)

  # Find neighboring genes on the left and right that have similar mutation rates
  left_genes <- get_genes_within_range(chromosome, leftBorder = leftmost_coordinate, rightBorder = left_targetGene_border) %>% base::intersect(genesSimMutRate)
  right_genes <- get_genes_within_range(chromosome, leftBorder = right_targetGene_border, rightBorder = rightmost_coordinate) %>% base::intersect(genesSimMutRate)

  # Combine left and right neighboring genes
  neighbor_genes <- c(left_genes, right_genes)

  # Handle the case when no neighboring genes are found
  if (base::length(neighbor_genes) == 0) {
    message("No neighboring genes satisfying the mutation rate similarity criteria were found within the specified window around ", hgnc, ". Expanding search to closest genes.")

    # Retrieve all genes with similar mutation rates from the GTF file
    gtf_genes <- mutmatch_gtf_grch37[mutmatch_gtf_grch37$type == "gene" &
      mutmatch_gtf_grch37$gene_id %in% genesSimMutRate &
      mutmatch_gtf_grch37$gene_id %in% mutmatch_topTranscripts$Gene]

    if (base::length(gtf_genes) > 0) {
      # Find closest neighboring genes on both sides if none are found initially
      if (base::length(left_genes) == 0) {
        left_genes <- find_closest_neighbor(grGene = grGene, grGenes = gtf_genes, side = "left")
      }
      if (base::length(right_genes) == 0) {
        right_genes <- find_closest_neighbor(grGene = grGene, grGenes = gtf_genes, side = "right")
      }
    }

    neighbor_genes <- c(left_genes, right_genes)
  }

  # If neighboring genes are found, get their genomic coordinates
  if (base::length(neighbor_genes) != 0) {
    coordinates <- base::lapply(neighbor_genes, function(x) {
      out <- get_gene_location(ensembl_gene_id = x)
      out[, gene := x]
      out[, hgnc := base::unique(mutmatch_topTranscripts[Gene == x]$Hugo_Symbol)]
      return(out)
    })
    coordinates <- data.table::rbindlist(coordinates)
    return(coordinates)
  } else {
    return(NULL)
  }
}
