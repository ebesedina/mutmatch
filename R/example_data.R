#' Example Mutations Data
#'
#' A compressed CSV file containing example somatic mutations.
#'
#' @format A compressed CSV file with columns:
#' \describe{
#'   \item{Sample}{Identifier for the sample where the mutation was observed.}
#'   \item{CHROM}{Chromosome where the mutation occurred.}
#'   \item{POS}{Genomic position of the mutation.}
#'   \item{REF}{Reference allele at the mutation site.}
#'   \item{ALT}{Mutated (alternate) allele at the mutation site.}
#' }
#' @name example_mutations.csv.gz
#' @docType data
NULL

#' Example Genome-wide Annotation Data
#'
#' A compressed CSV file containing example genome-wide annotation information.
#' It should have at least two columns: Sample and Feature.
#'
#' @format A compressed CSV file with columns:
#' \describe{
#'   \item{Sample}{Identifier for the sample.}
#'   \item{Feature1}{Sample annotation, e.g., cancer type, treatment, etc.}
#'   \item{Feature2}{Sample annotation, e.g., cancer type, treatment, etc.}
#' }
#' @name example_genomewide_annotation.csv.gz
#' @docType data
NULL

#' Example Gene Annotation Data
#'
#' A compressed CSV file containing example gene annotation information.
#'
#' @format A compressed CSV file with columns:
#' \describe{
#'   \item{Sample}{Identifier for the sample.}
#'   \item{Gene}{Start position of the gene.}
#'   \item{Feature}{Gene feature, for example copy number state.}
#' }
#' @name example_gene_annotation.csv.gz
#' @docType data
NULL
