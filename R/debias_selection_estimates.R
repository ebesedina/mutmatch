#' Debias Selection Estimates
#'
#' This function debiases selection estimates by taking the original estimates
#' and comparing them with simulated estimates. It corrects the coefficients
#' based on the median value of the simulated coefficients.
#'
#' @param estimates_original A data.table containing the original coefficient estimates.
#' @param estimates_simulated A data.table containing the simulated coefficient estimates.
#'
#' @return A data.table with corrected coefficient estimates for the isTarget variable and associated terms.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' estimates_original <- data.table::data.table(
#'   coefName = c("isTarget", "CNA"),
#'   coef = c(0.1, 0.2)
#' )
#' estimates_simulated <- data.table::data.table(coefName = rep(c("isTarget", "CNA"),
#'   each = 4
#' ), coef = runif(8))
#' debias_selection_estimates(estimates_original, estimates_simulated)
#' }
debias_selection_estimates <- function(estimates_original, estimates_simulated) {
  # Add an index column to remember the original order of estimates_original
  estimates_original[, orig_order := seq_len(.N)]

  if (!"Cohort" %in% colnames(estimates_original)) {
    # Summarize the distribution of simulated coefficients by taking the median
    estimates_simulated <- estimates_simulated %>%
      dplyr::group_by(coefName) %>%
      dplyr::summarise(coef_sim_median = median(coef)) %>%
      data.table::data.table()

    # Merge the original and simulated estimates based on coefficient names
    estimates_corrected <- merge(estimates_original, estimates_simulated, by = c("coefName"))
  } else {
    # Summarize the distribution of simulated coefficients by taking the median
    estimates_simulated <- estimates_simulated %>%
      dplyr::group_by(coefName, Cohort) %>%
      dplyr::summarise(coef_sim_median = median(coef)) %>%
      data.table::data.table()

    # Merge the original and simulated estimates based on coefficient names
    estimates_corrected <- merge(estimates_original, estimates_simulated, by = c("coefName", "Cohort"))
  }

  # Sort by the original order
  data.table::setorder(estimates_corrected, orig_order)

  # Remove the temporary index column
  estimates_corrected[, orig_order := NULL]

  # Correct the coefficients
  estimates_corrected[, coef_corrected := ifelse(stringr::str_detect(coefName, "isTarget"),
    coef - coef_sim_median,
    coef
  )]

  # Reorder columns to make sure coef_corrected is the third column
  data.table::setcolorder(estimates_corrected, c(
    "coefName", "coef", "coef_corrected",
    setdiff(
      names(estimates_corrected),
      c("coefName", "coef", "coef_corrected")
    )
  ))

  return(estimates_corrected[])
}
