#' Get CADD Regression Table
#'
#' Creates a mutation table for a specified target gene, categorizing mutations based on their
#' Combined Annotation Dependent Depletion (CADD) scores into high-CADD and low-CADD regions.
#' It also performs filtering based on regions.
#'
#' @inheritParams get_gene_neighbors
#' @inheritParams filter_ranges
#' @inheritParams load_annotation
#' @inheritParams load_mutations
#' @param caddScoresPath Specifies the file path to a bigWig-formatted file containing CADD scores for the hg19 genome.
#' The default version (CADD_GRCh37-v1.4.bw) can be downloaded using the function \code{\link{download_cadd_file}}.
#'
#' @return A data.table containing categorized mutation information for the target gene and its neighboring genes.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' mutation_table <- get_cadd_model_table(hgnc = "BRCA1",
#' caddScoresPath = "path/to/CADD_scores.bw", ...)
#' }
#'
#' @export
get_cadd_model_table <- function(hgnc,
                                 caddScoresPath,
                                 mutationsPath,
                                 annotationGenePath,
                                 annotationGenomeWidePath,
                                 filterRegionsPath = NULL,
                                 format = NA,
                                 filter = NA) {
  # Retrieve the Ensembl ID corresponding to the given HGNC symbol
  ensembl_id <- get_ensembl_from_hgnc(hgnc)

  # Get the gene's genomic location data and create genomic range objects
  geneData <- get_gene_location(ensembl_gene_id = ensembl_id)
  grGene <- make_ranges_from_table(geneData, exonsCoordinates = FALSE)
  grTargetGene <- make_ranges_from_table(geneData, exonsCoordinates = TRUE)

  # Import CADD data
  CADD_local <- rtracklayer::import(con = caddScoresPath, which = grTargetGene) %>% suppressWarnings()

  # Define CADD score threshold
  threshold_constr <- 20

  # Pre-processing to avoid conflicts in chromosome length
  GenomeInfoDb::seqlengths(CADD_local) <- NA

  # Split data based on CADD score threshold
  grTargetGene_constrained <- CADD_local[CADD_local$score >= threshold_constr]
  grTargetGene_unconstrained <- CADD_local[CADD_local$score < threshold_constr]

  # Intersect the constrained and unconstrained regions with the target gene
  grTargetGene_constrained <- GenomicRanges::intersect(grTargetGene, grTargetGene_constrained, ignore.strand = TRUE) %>%
    GenomicRanges::reduce()
  grTargetGene_unconstrained <- GenomicRanges::intersect(grTargetGene, grTargetGene_unconstrained, ignore.strand = TRUE) %>%
    GenomicRanges::reduce()

  # Filter unwanted regions from the genomic ranges of the target gene parts
  grTargetGene_constrained <- filter_ranges(
    grObject = grTargetGene_constrained,
    filterRegionsPath = filterRegionsPath,
    format = format,
    filter = filter
  )
  grTargetGene_unconstrained <- filter_ranges(
    grObject = grTargetGene_unconstrained,
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

  # Generate mutation tables for the target gene low- and high-CADD regions
  targetMutationTable <- get_mutations_table(
    grObject = grTargetGene_constrained,
    grMutations = grMutations,
    annotation = annotation
  )
  baselineMutationTable <- get_mutations_table(
    grObject = grTargetGene_unconstrained,
    grMutations = grMutations,
    annotation = annotation
  )

  # Mark rows to differentiate high-CADD and low-CADD regions
  targetMutationTable$isTarget <- factor(1)
  baselineMutationTable$isTarget <- factor(0)

  # Combine mutation tables
  mutation_table <- rbind(targetMutationTable, baselineMutationTable)

  # Relevel factors in the isTarget column
  mutation_table$isTarget <- stats::relevel(mutation_table$isTarget, ref = "0")

  return(mutation_table)
}
