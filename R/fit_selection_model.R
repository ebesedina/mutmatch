#' Fit Selection Model and Debias Its Estimates
#'
#' This function fits a selection model to provided mutation data based on a user-defined formula.
#' The process involves an initial simulation of the data under a neutral selection model.
#' Following this, the chosen model is applied to both the real and simulated datasets.
#' The output estimates are then corrected using the simulated data as a reference.
#' Additionally, the function aggregates and provides the total number of mutations
#' for both tested and controlled regions across distinct cohorts.
#'
#' @inheritParams fit_glm
#' @param mutation_table A data.table containing the mutation data. The table should
#'   include mandatory columns 'MutationNumber' and 'isTarget', as well as any additional
#'   covariates that may be specified in the model formula. If any covariates are present
#'   in the data but not specified in the formula, the data will be subsetted according to
#'   these covariates to fit individual models for each subset.
#' @param simtimes The number of times for bias correction. Defaults to 50.
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
    message("Performing model fitting on real data")
    estimates_original <- fit_glm(data = mutation_table, formula = formula, family = family)

    # Fit the model on the simulated data
    message("Performing model fitting on simulated data")
    estimates_simulated <- pbapply::pblapply(1:simtimes, function(i) {
      fit_glm(data = mutation_table_shuffled[permutation_id == i, ], formula = formula)
    }) %>% data.table::rbindlist()

    # Correct the estimates
    message("Correcting selection estimates using simulated data")
    estimates_corrected <- debias_selection_estimates(estimates_original, estimates_simulated)

    # Add mutation counts per tested and controlled region
    mut_counts_total <- mutation_table[, .(MutationNumber = sum(MutationNumber)), by = .(isTarget)]
    mut_counts_total <- data.table::dcast(mut_counts_total, formula = ... ~ isTarget, value.var = "MutationNumber")
    data.table::setnames(mut_counts_total, old = c("0", "1"), new = c("MutationNumber_control", "MutationNumber_tested"))

    # Remove the column named "."
    mut_counts_total[, "." := NULL]

    estimates_corrected <- cbind(estimates_corrected, mut_counts_total)
  } else {
    # If there are model-specific cohorts
    data.table::setDT(mutation_table)[, Cohort := do.call(paste, c(.SD, sep = "__")), .SDcols = model_cohorts_vars]
    data.table::setDT(mutation_table_shuffled)[, Cohort := do.call(paste, c(.SD, sep = "__")), .SDcols = model_cohorts_vars]
    model_cohorts <- unique(mutation_table$Cohort)

    # Fit the model on the real data for each cohort
    message("Performing model fitting on real data for each cohort")
    estimates_original <- pbapply::pblapply(model_cohorts, function(cohort) {
      mutation_table_cohort <- mutation_table[Cohort == cohort]
      estimates_original_cohort <- fit_glm(data = mutation_table_cohort, formula = formula)
      estimates_original_cohort[, Cohort := cohort]
      estimates_original_cohort
    }) %>% data.table::rbindlist()

    # Fit the model on the simulated data for each cohort
    message("Performing model fitting on simulated data for each cohort")
    # Create a data frame of all combinations of cohorts and simtimes
    combinations <- data.table::CJ(Cohort = model_cohorts, Simtime = 1:simtimes)

    # Loop through each row of combinations using pblapply
    estimates_simulated <- pbapply::pblapply(seq_len(nrow(combinations)), function(index) {
      cohort <- combinations$Cohort[index]
      i <- combinations$Simtime[index]
      mutation_table_shuffled_i <- mutation_table_shuffled[Cohort == cohort &
        permutation_id == i, ]
      estimates_simulated_cohort <- fit_glm(data = mutation_table_shuffled_i, formula = formula)
      estimates_simulated_cohort[, Cohort := cohort]
      return(estimates_simulated_cohort)
    }) %>% data.table::rbindlist()

    # Correct the estimates
    message("Correcting selection estimates using simulated data")
    estimates_corrected <- debias_selection_estimates(estimates_original, estimates_simulated)

    # Add mutation counts per tested and controlled region
    mut_counts_total <- mutation_table[, .(MutationNumber = sum(MutationNumber)), by = .(isTarget, Cohort)]
    mut_counts_total <- data.table::dcast(mut_counts_total, formula = Cohort ~ isTarget, value.var = "MutationNumber")
    data.table::setnames(mut_counts_total, old = c("0", "1"), new = c("MutationNumber_control", "MutationNumber_tested"))
    estimates_corrected <- merge(estimates_corrected, mut_counts_total, by = "Cohort")
  }

  return(estimates_corrected)
}
