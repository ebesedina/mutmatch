#' Remove Mismatched Samples
#'
#' This function removes samples from a reference annotation table whose gene-specific annotations
#' do not match those in the provided list of genes.
#'
#' @inheritParams load_annotation
#' @param referenceAnnotation Reference annotation from which mismatched samples have to be removed.
#' @param genes A vector of HGNC gene symbols to be used for filtering samples.
#'
#' @return A modified version of the `referenceAnnotation` data.table or data.frame, containing only samples
#'         that match the gene-specific annotations.
#'
#' @examples
#' # Usage:
#' # filtered_reference = remove_mismatched_samples(genes, annotationGenePath, referenceAnnotation)
#' @export
remove_mismatched_samples <- function(genes, annotationGenePath, referenceAnnotation) {
  # Initialize an empty list to store matched samples for each gene
  matched_neighbors_samples <- lapply(genes, function(hgnc_i) {
    # Load gene-specific annotation
    annotation_gene <- load_annotation(
      annotationGenePath = annotationGenePath,
      annotationGenomeWidePath = NULL,
      hgnc = hgnc_i
    )

    # Check if gene-specific annotation exists
    if (nrow(annotation_gene) > 0) {
      # Merge with reference annotation and extract unique sample names
      annotation_matched <- base::merge(annotation_gene, referenceAnnotation, by = colnames(annotation_gene))
      matched_samples <- unique(annotation_matched$Sample)
    } else {
      # If gene has no gene-specific annotation, assume it has the same annotation as reference
      matched_samples <- unique(referenceAnnotation$Sample)
    }
    return(matched_samples)
  })

  # Intersect all matched samples lists to find common samples
  matched_samples <- base::Reduce(base::intersect, matched_neighbors_samples)

  # Subset reference annotation to include only matched samples
  referenceAnnotation <- referenceAnnotation[Sample %in% matched_samples]

  return(referenceAnnotation)
}
