#' Fit a selection model and debias its estimates
#'
#' This function fits a selection model to given mutation data using a specified formula.
#' The function first simulates the data using a neutral selection model and then fits
#' the specified model to both the original and simulated data. It corrects the estimates
#' for selection using the simulated data as a reference.
#'
#' @inheritParams fit_glm
#' @param mutation_table A data.table containing the mutation data. The table should
#'   include mandatory columns 'MutationNumber' and 'isTarget', as well as any additional
#'   covariates that may be specified in the model formula. If any covariates are present
#'   in the data but not specified in the formula, the data will be subsetted according to
#'   these covariates to fit individual models for each subset.
#' @param simtimes The number of times to simulate the data for correction of selection
#'   estimates. Defaults to 50.
#'
#' @return A data.table containing corrected estimates.
#' @export
fit_selection_model <- function(mutation_table,
                                formula = "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))",
                                family = "bayes.poisson",
                                simtimes = 50) {
  # Message for tracking progress
  message("Simulating data using neutral selection model")

  # Identify the cohort variables in the data
  cohort_vars <- setdiff(names(mutation_table), c("isTarget", "Mutation", "MutationNumber", "ntAtRisk"))

  # Shuffle mutation counts separately for each cohort
  data.table::setDT(mutation_table)[, Cohort := do.call(paste, c(.SD, sep = "__")), .SDcols = cohort_vars]
  cohort_ids <- unique(mutation_table$Cohort)
  mutation_table_shuffled <- lapply(cohort_ids, function(cohort_i) {
    shuffle_mutation_counts(mutation_table[Cohort == cohort_i, -"Cohort", with = FALSE], ntimes = simtimes)
  }) %>% data.table::rbindlist()

  # Identify the model-specific cohorts
  mutation_table$Cohort <- NULL
  model_cohorts_vars <- setdiff(names(mutation_table), extract_variables_from_formula(formula))

  # If there are no model-specific cohorts
  if (length(model_cohorts_vars) == 0) {
    # Fit the model on the real data
    message("Fitting the model in the real data")
    estimates_original <- fit_glm(data = mutation_table, formula = formula, family = family)

    # Fit the model on the simulated data
    message("Fitting the model in the simulated data")
    estimates_simulated <- pbapply::pblapply(1:simtimes, function(i) {
      fit_glm(data = mutation_table_shuffled[permutation_id == i, ], formula = formula)
    }) %>% data.table::rbindlist()

    # Correct the estimates
    message("Correcting selection estimates using simulated data")
    estimates_corrected <- debias_selection_estimates(estimates_original, estimates_simulated)
  } else {
    # If there are model-specific cohorts
    data.table::setDT(mutation_table)[, Cohort := do.call(paste, c(.SD, sep = "__")), .SDcols = model_cohorts_vars]
    data.table::setDT(mutation_table_shuffled)[, Cohort := do.call(paste, c(.SD, sep = "__")), .SDcols = model_cohorts_vars]
    model_cohorts <- unique(mutation_table$Cohort)

    # Fit the model on the real data for each cohort
    message("Fitting the model in the real data")
    estimates_original <- lapply(model_cohorts, function(cohort) {
      mutation_table_cohort <- mutation_table[Cohort == cohort]
      estimates_original_cohort <- fit_glm(data = mutation_table_cohort, formula = formula)
      estimates_original_cohort[, Cohort := cohort]
      estimates_original_cohort
    }) %>% data.table::rbindlist()

    # Fit the model on the simulated data for each cohort
    message("Fitting the model in the simulated data")
    estimates_simulated <- lapply(model_cohorts, function(cohort) {
      pbapply::pblapply(1:simtimes, function(i) {
        mutation_table_shuffled_i <- mutation_table_shuffled[Cohort == cohort & permutation_id == i, ]
        estimates_simulated_cohort <- fit_glm(data = mutation_table_shuffled_i, formula = formula)
        estimates_simulated_cohort[, Cohort := cohort]
        estimates_simulated_cohort
      }) %>% data.table::rbindlist()
    }) %>% data.table::rbindlist()

    # Correct the estimates
    message("Correcting selection estimates using simulated data")
    estimates_corrected <- debias_selection_estimates(estimates_original, estimates_simulated)
  }

  return(estimates_corrected)
}
