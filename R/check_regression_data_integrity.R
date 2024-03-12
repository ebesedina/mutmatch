#' Check Regression Data Integrity
#'
#' This function checks the integrity of the data to be used in a regression model. It ensures that
#' the `isTarget` variable has two levels and checks the presence of all combinations of `isTarget`
#' and an interaction variable from the formula in the dataset. If any combination is missing, it
#' will throw an error.
#'
#' @param mutation_table A data.table containing the dataset to check.
#' @param formula A character string representing the formula used in the regression model.
#'        Default is "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))".
#'
#' @return This function does not return a value. It will stop and throw an error message
#'         if the data integrity checks fail.
#' @examples
#' \dontrun{
#' check_regression_data_integrity(mutation_table)
#' }
#' @export
check_regression_data_integrity <- function(mutation_table,
                                            formula = "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))") {
  # Check presence of columns used in the formula
  var_names <- extract_variables_from_formula(formula)
  non_existing_vars <- setdiff(var_names, names(mutation_table))
  if (length(non_existing_vars) > 0) {
    stop(sprintf(
      "The following variables from the formula are not found in the data.table: %s",
      paste(non_existing_vars, collapse = ", ")
    ))
  }

  # Check if isTarget has 2 levels
  if (length(unique(mutation_table$isTarget)) < 2) {
    if (all(mutation_table$isTarget == 0)) {
      stop("There is no DNA in the tested region for isTarget variable, not able to fit the model")
    } else {
      stop("There is no DNA in the control region for isTarget variable, not able to fit the model")
    }
  }
  # Extract interaction variable from the formula
  interaction_var <- gsub(".*isTarget:([a-zA-Z0-9_]+).*", "\\1", formula)

  if (interaction_var == formula) interaction_var <- NULL # No interaction variable

  if (!is.null(interaction_var)) {
    if (length(unique(mutation_table[[interaction_var]])) < 2) {
      if (all(mutation_table[[interaction_var]]) == 0) {
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
    present_combinations <- unique(mutation_table[, .(isTarget, get(interaction_var))])

    data.table::setnames(all_combinations, "V2", interaction_var) # Rename the column
    data.table::setnames(present_combinations, "V2", interaction_var) # Rename the column

    missing_combinations <- all_combinations[!present_combinations, on = c("isTarget", interaction_var)]

    if (nrow(missing_combinations) > 0) {
      missing_combinations_str <- apply(missing_combinations, 1, function(x) {
        paste(names(missing_combinations), x, sep = ": ", collapse = ", ")
      })
      missing_combinations_str <- paste(missing_combinations_str, collapse = "; ")
      stop(sprintf("Missing combinations for isTarget:%s: %s", interaction_var, missing_combinations_str))
    }
  }
}
