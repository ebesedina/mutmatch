#' Get Neighbors Regression Table
#'
#' This function generates a mutation table for a target gene and its neighboring genes.
#' It also performs filtering based on regions and annotations.
#'
#' @inheritParams get_gene_neighbors
#' @inheritParams filter_ranges
#' @inheritParams load_annotation
#' @inheritParams load_mutations
#'
#' @return A data.table containing mutation information for the target gene and its neighbors.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' mutation_table <- get_neighbors_model_table(...)
#' }
#' @export
get_neighbors_model_table <- function(hgnc,
                                      mutationsPath,
                                      annotationGenePath,
                                      annotationGenomeWidePath,
                                      filterRegionsPath = NULL,
                                      format = NA,
                                      filter = NA,
                                      neighborsWindow = "0.5Mb",
                                      outlierNeighborsThreshold = 0.2) {
  # Retrieve the Ensembl ID corresponding to the given HGNC symbol
  ensembl_id <- get_ensembl_from_hgnc(hgnc)

  # Get the gene's genomic location data and create genomic range objects
  geneData <- get_gene_location(ensembl_gene_id = ensembl_id)
  grGene <- make_ranges_from_table(geneData, exonsCoordinates = FALSE)
  grTargetGene <- make_ranges_from_table(geneData, exonsCoordinates = TRUE)

  # Retrieve neighboring genes within a given window around the target gene
  neighborsData <- get_gene_neighbors(
    grGene = grGene,
    hgnc = hgnc,
    neighborsWindow = neighborsWindow,
    outlierNeighborsThreshold = outlierNeighborsThreshold
  )

  if (is.null(neighborsData)) stop(base::paste0("There are no neighboring genes for ", hgnc, " with given parameters"))

  # Create genomic ranges for neighboring genes
  grNeighbours <- make_ranges_from_table(neighborsData, exonsCoordinates = TRUE)

  # Filter unwanted regions from the genomic ranges of the target gene and neighbors
  grTargetGene <- filter_ranges(
    grObject = grTargetGene,
    filterRegionsPath = filterRegionsPath,
    format = format,
    filter = filter
  )

  # Code for filtering ranges of neighboring genes (only CUP and CRG75)
  grNeighbours <- filter_ranges(
    grObject = grNeighbours,
    filterRegionsPath = NULL,
    format = NA,
    filter = NA
  )

  # Load mutations and annotation data for the chromosome of interest
  chr <- GenomeInfoDb::seqnames(grGene) %>% base::as.character()
  grMutations <- load_mutations(mutationsPath = mutationsPath, chr = chr)
  annotation <- load_annotation(
    annotationGenePath = annotationGenePath,
    annotationGenomeWidePath = annotationGenomeWidePath,
    hgnc = hgnc
  )

  # Remove samples that don't match annotation between target and neighboring genes
  if (!is.null(annotationGenePath)) {
    neighbors_hgnc <- neighborsData$hgnc %>% base::unique()
    annotation <- remove_mismatched_samples(
      genes = neighbors_hgnc,
      annotationGenePath = annotationGenePath,
      referenceAnnotation = annotation
    )
  }

  # Generate mutation tables for the target gene and its neighbors
  targetGeneMutationTable <- get_mutations_table(
    grObject = grTargetGene,
    grMutations = grMutations,
    annotation = annotation
  )

  neighborsMutationTable <- get_mutations_table(
    grObject = grNeighbours,
    grMutations = grMutations,
    annotation = annotation
  )

  # Mark rows to differentiate target gene from its neighbors
  targetGeneMutationTable$isTarget <- factor(1)
  neighborsMutationTable$isTarget <- factor(0)

  # Combine mutation tables
  mutation_table <- rbind(targetGeneMutationTable, neighborsMutationTable)

  # Relevel factors in the isTarget column
  mutation_table$isTarget <- relevel(mutation_table$isTarget, ref = "0")
  return(mutation_table)
}
