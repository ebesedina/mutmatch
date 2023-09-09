#' Extract Variables from a Formula String
#'
#' This function takes a formula string as input and returns a character vector containing all the variables in the formula.
#' The function does not remove any offset terms if present.
#'
#' @param formula A character string representing the formula.
#'
#' @return A character vector containing all the variables in the formula.
#'
#' @examples
#' \dontrun{
#' formula_str <- "MutationNumber ~ isTarget + CNA + Mutation + offset(log(ntAtRisk))"
#' extracted_vars <- extract_variables_from_formula(formula_str)
#' print(extracted_vars)
#' }
#'
#' @export
extract_variables_from_formula <- function(formula) {
  # Convert the formula string to a formula object
  formula <- stats::as.formula(formula)

  # Extract all variable names
  all_vars <- base::all.vars(formula)

  return(all_vars)
}
