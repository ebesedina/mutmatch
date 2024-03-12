#' Correct P-Values to Remove Deflation or Inflation
#'
#' @description This function corrects for deflation or inflation in real p-values by mapping them to their
#' corresponding quantiles within a distribution of simulated p-values. This adjustment is based on the empirical
#' cumulative distribution function (ECDF) derived from the simulated p-values, enhancing the statistical reliability
#' of the p-value estimates.
#' @param real_pvalues A numeric vector of real p-values to be adjusted.
#' @param simulated_pvalues A numeric vector of simulated p-values used to generate the ECDF.
#' @return A numeric vector of adjusted p-values.
#' @examples
#' \dontrun{
#' real_pvalues <- c(0.01, 0.05, 0.2)
#' simulated_pvalues <- runif(1000, 0, 1)
#' corrected_pvalues <- correct_pvalues(real_pvalues, simulated_pvalues)
#' }
#' @export
correct_pvalues <- function(real_pvalues, simulated_pvalues) {
  ecdf_simulated <- ecdf(simulated_pvalues)
  adjusted_pvalues <- ecdf_simulated(real_pvalues)
  return(adjusted_pvalues)
}

#' Process Target Coefficient P-Values
#'
#' @description For specific coefficients identified by `target_type`, this function adjusts their p-values using
#' the `correct_pvalues` function. It is applied to a subset of data where the coefficient names match `target_type`.
#' @param data A data frame containing the selection estimates, including both real and simulated p-values.
#' @param target_type A character string specifying the coefficient name to filter and process the p-values for.
#' @return The data frame with corrected p-values for the specified coefficient.
#' @examples
#' \dontrun{
#' processed_data <- process_target_rows(selection_estimates, "isTarget1")
#' }
#' @export
process_target_rows <- function(data, target_type) {
  target_rows <- data$coefName == target_type
  if (any(target_rows)) {
    pvalues_strings <- as.character(data[target_rows, "pVal_neglog10_sim"])

    # Use stringr::str_extract_all to find all matches of the pattern (numeric values including decimals)
    # This returns a list of character vectors
    extracted_pvalues_list <- stringr::str_extract_all(pvalues_strings, "\\d+(?:\\.\\d+)?")
    pvalues_sim <- 10^(-as.numeric(unlist(extracted_pvalues_list)))
    data[target_rows, pVal_corrected := correct_pvalues(data[target_rows]$pVal, pvalues_sim)]
  }
  return(data[])
}


#' Correct P-Values for Selection Estimates Across Multiple Genes
#'
#' @description This function refines p-value estimates for selection analyses by correcting biases
#' across a comprehensive range of genes. It applies the `correct_pvalues` function to adjust p-values for
#' specific coefficients (i.e., `isTarget1`, `isTarget1:CNA-1`, and `isTarget1:CNA1`) in the selection estimates dataset.
#' For enhanced accuracy, it's recommended to apply this function to data encompassing a broad spectrum of genes.
#' Alternatively, increase the simulation count (`simtimes` parameter) in functions like `get_selection_estimates_neighbors`
#' or `get_selection_estimates_cadd` to improve resolution when extensive gene data is unavailable.
#' @param selection_estimates A data frame of selection estimates that includes columns for real and simulated p-values.
#' @return A data frame with corrected p-values for the specified coefficients, enhancing the reliability of selection analyses.
#' @examples
#' \dontrun{
#' selection_estimates <- get_selection_estimates_neighbors(
#'   hgnc = "KRAS",
#'   mutationsPath = system.file("extdata", "example_mutations.csv.gz", package = "mutmatch"),
#'   annotationGenePath = system.file("extdata", "example_gene_annotation.csv.gz",
#'     package = "mutmatch"
#'   ),
#'   annotationGenomeWidePath = system.file("extdata", "example_genomewide_annotation.csv",
#'     package = "mutmatch"
#'   ),
#'   neighborsWindow = "0.5Mb",
#'   outlierNeighborsThreshold = 0.2,
#'   simtimes = 500
#' )
#' selection_estimates <- post_process_pvalues(selection_estimates)
#' }
#' @export
post_process_pvalues <- function(selection_estimates) {
  target_coefficients <- c("isTarget1", "isTarget1:CNA-1", "isTarget1:CNA1")
  if ("pVal_neglog10_sim" %in% colnames(selection_estimates)) {
    for (target in target_coefficients) {
      selection_estimates <- process_target_rows(selection_estimates, target)
    }
  }
  return(selection_estimates[])
}
