check_regression_data_integrity <- function(mutation_table,
                                            formula = "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))") {
  # Check if isTarget has 2 levels
  if (length(unique(mutation_table$isTarget)) < 2) {
    if (unique(mutation_table$isTarget) == 0) {
      stop("There is no DNA in the tested region for isTarget variable, not able to fit the model")
    } else {
      stop("There is no DNA in the control region for isTarget variable, not able to fit the model")
    }
  }
  # Extract interaction variable from the formula
  interaction_var <- gsub(".*isTarget:([a-zA-Z_]+).*", "\\1", formula)
  if (interaction_var == formula) interaction_var <- NULL # No interaction variable

  if (!is.null(interaction_var)) {
    if (length(unique(mutation_table[[interaction_var]])) < 2) {
      if (unique(mutation_table[[interaction_var]]) == 0) {
        stop(paste0("There is no DNA in the tested region for", interaction_var, "variable, not able to fit the model"))
      } else {
        stop(paste0("There is no DNA in the control region for", interaction_var, "variable, not able to fit the model"))
      }
    }
    # Check if all combinations of isTarget and interaction_var are present
    all_combinations <- data.table::CJ(
      isTarget = c("0", "1"),
      V2 = unique(mutation_table[[interaction_var]])
    )
    data.table::setnames(all_combinations, "V2", interaction_var) # Rename the column

    present_combinations <- unique(mutation_table[, .(isTarget, get(interaction_var))])
    data.table::setnames(present_combinations, "V2", interaction_var) # Rename the column

    if (nrow(all_combinations) != nrow(present_combinations)) {
      missing_combinations <- all_combinations[!present_combinations, on = .(isTarget, get(interaction_var))]
      stop(sprintf("Missing combinations for isTarget:%s: %s", interaction_var, toString(missing_combinations)))
    }
  }
}
