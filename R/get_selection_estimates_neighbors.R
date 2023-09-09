#' Estimate Selection Forces on Neighboring Genes
#'
#' This function performs a two-step analysis to estimate the selection pressures acting on a given target gene and its neighboring genes.
#' The function first identifies neighboring genes within a specified genomic window and filters them based on a mutation rate similarity threshold.
#' Then, it fits a Bayesian Poisson regression model to calculate the debiased selection estimates.
#'
#' @inheritParams get_neighbors_model_table
#' @inheritParams fit_selection_model
#'
#' @return A data.table containing the debiased selection estimates, including estimated coefficients and their uncertainties.
#'
#' @examples
#' \dontrun{
#' selection_estimates <- get_selection_estimates_neighbors(
#'   hgnc = "KRAS",
#'   mutationsPath = system.file("extdata",
#'     "example_mutations.csv.gz",
#'     package = "mutmatch"
#'   ),
#'   annotationGenePath = system.file("extdata", "example_gene_annotation.csv.gz",
#'     package = "mutmatch"
#'   ),
#'   annotationGenomeWidePath = system.file("extdata",
#'     "example_genomewide_annotation.csv",
#'     package = "mutmatch"
#'   ),
#'   neighborsWindow = "0.5Mb",
#'   outlierNeighborsThreshold = 0.2
#' )
#' }
#' @export
get_selection_estimates_neighbors <- function(hgnc,
                                              mutationsPath,
                                              annotationGenePath,
                                              annotationGenomeWidePath,
                                              filterRegionsPath = NULL,
                                              format = NA,
                                              filter = NA,
                                              neighborsWindow = "0.5Mb",
                                              outlierNeighborsThreshold = 0.2,
                                              formula = "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))",
                                              family = "bayes.poisson",
                                              simtimes = 50) {
  # Get the data table for the regression model
  message("Retrieving mutation counts using neighbors baseline mutation model")
  mutation_table <- get_neighbors_model_table(
    hgnc = hgnc,
    mutationsPath = mutationsPath,
    annotationGenePath = annotationGenePath,
    annotationGenomeWidePath = annotationGenomeWidePath,
    filterRegionsPath = filterRegionsPath,
    format = format,
    filter = filter,
    neighborsWindow = neighborsWindow,
    outlierNeighborsThreshold = outlierNeighborsThreshold
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
