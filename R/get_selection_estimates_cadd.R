#' Estimate Selection Pressures on a Target Gene Using Low-CADD Gene Regions as a Mutational Rate Baseline
#'
#' Performs a comprehensive two-stage analysis to quantify the selection forces acting upon a specified target gene.
#' In the first stage, the function separates gene regions with low-CADD and high-CADD scores.
#' In the second stage, a Bayesian Poisson regression model is employed to derive selection estimates for the target gene.
#' These estimates are subsequently refined through permutation-based debiasing.

#' @inheritParams get_cadd_model_table
#' @inheritParams fit_selection_model
#'
#' @return A data.table containing the debiased selection estimates, including estimated coefficients and their uncertainties.
#'
#' @examples
#' \dontrun{
#' # Replace "your/destination/path/CADD_GRCh37-v1.4.bw" with the path where
#' you want to save the file.
#'
#' caddScoresPath = "your/destination/path/CADD_GRCh37-v1.4.bw"
#' download_cadd_file(caddScoresPath = caddScoresPath)
#'
#' selection_estimates <- get_selection_estimates_cadd(
#'   hgnc = "KRAS",
#'   caddScoresPath = caddScoresPath,
#'   mutationsPath = system.file("extdata",
#'     "example_mutations.csv.gz",
#'     package = "mutmatch"
#'   ),
#'   annotationGenePath = system.file("extdata",
#'   "example_gene_annotation.csv.gz",
#'     package = "mutmatch"
#'   ),
#'   annotationGenomeWidePath = system.file("extdata",
#'     "example_genomewide_annotation.csv",
#'     package = "mutmatch"
#'   )
#' )
#' }
#' @export
get_selection_estimates_cadd <- function(hgnc,
                                         caddScoresPath,
                                         mutationsPath,
                                         annotationGenePath,
                                         annotationGenomeWidePath,
                                         filterRegionsPath = NULL,
                                         format = NA,
                                         filter = NA,
                                         formula = "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))",
                                         family = "bayes.poisson",
                                         simtimes = 50) {
  # Get the data table for the regression model
  message("Retrieving mutation counts using CADD baseline mutation model")
  mutation_table <- get_cadd_model_table(
    hgnc = hgnc,
    caddScoresPath = caddScoresPath,
    mutationsPath = mutationsPath,
    annotationGenePath = annotationGenePath,
    annotationGenomeWidePath = annotationGenomeWidePath,
    filterRegionsPath = filterRegionsPath,
    format = format,
    filter = filter
  )

  # Fit the selection model and get the debiased estimates
  message("Fitting the selection model and debiasing estimates")
  selection_estimates <- fit_selection_model(mutation_table,
    formula = formula,
    family = family,
    simtimes = simtimes
  )

  return(selection_estimates)
}
