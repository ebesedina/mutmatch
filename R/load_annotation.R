#' Load the file(s) with annotation
#'
#' @param annotationGenomeWidePath A path to the file with genome-wide annotation for each sample
#' @param annotationGenePath A path to the file with gene-specific annotation for each sample. Default is NULL.
#' @param hgnc The HGNC symbol of the target gene in case a gene-specific annotation is used. Default is NULL.
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \dontrun{
#' annotationGenomeWidePath <- system.file("extdata",
#'   "example_genomewide_annotation.csv",
#'   package = "mutmatch"
#' )
#' annotationGenePath <- system.file("extdata", "example_gene_annotation.csv.gz",
#'   package = "mutmatch"
#' )
#' hgnc <- "KRAS"
#' }
load_annotation <- function(annotationGenomeWidePath,
                            annotationGenePath = NULL,
                            hgnc = NULL) {
  if (!base::is.null(annotationGenomeWidePath)) {
    sample_annotation_genome <- data.table::fread(annotationGenomeWidePath) %>% unique()
  }
  if (!base::is.null(annotationGenePath)) {
    # Retrieve possible secondary names for the gene of interest
    gene_names <- check_gene_aliases(hgnc = hgnc, approved = FALSE)

    # Get only rows with gene name in it
    sample_annotation_gene <- data.table::fread(cmd = base::paste(
      "zcat",
      annotationGenePath,
      "| grep -E",
      base::paste0("'", base::paste(gene_names, collapse = "|"), "'")
    ))
    sample_annotation_gene_colnames <- base::scan(annotationGenePath, what = "", nlines = 1, sep = ",", quiet = T)
    colnames(sample_annotation_gene) <- sample_annotation_gene_colnames
    sample_annotation_gene <- sample_annotation_gene[Gene %in% gene_names][, -"Gene"] %>% unique()
  }

  # Return or both annotations or only the one which is specified
  if (!base::is.null(annotationGenomeWidePath) &
    !base::is.null(annotationGenePath)) {
    sample_annotation <- base::merge(sample_annotation_genome, sample_annotation_gene, by = "Sample") %>% unique()
  } else if (!base::is.null(annotationGenomeWidePath) &
    base::is.null(annotationGenePath)) {
    sample_annotation <- sample_annotation_genome %>% unique()
  } else if (base::is.null(annotationGenomeWidePath) &
    !base::is.null(annotationGenePath)) {
    sample_annotation <- sample_annotation_gene %>% unique()
  }

  return(sample_annotation)
}
