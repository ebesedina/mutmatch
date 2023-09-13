#' A Wrapper Function for fastglm::fastglmPure
#'
#' This function serves as a wrapper around fastglm::fastglmPure,
#' providing additional functionality and handling various input configurations.
#'
#' @importFrom fastglm fastglmPure
#'
#' @inheritParams fastglm::fastglmPure
#' @param control A list of control parameters including `maxit` for the maximum number of iterations (default is 100) and `epsilon` for the tolerance level (default is 1e-07).
#' @return A fitted generalized linear model.
#'
.fastglm_fit <- function(x,
                         y,
                         weights = rep(1, NROW(y)),
                         family = stats::gaussian(),
                         offset = rep(0, NROW(y)),
                         start = NULL,
                         mustart = NULL,
                         etastart = NULL,
                         control = list(maxit = 100, epsilon = 1e-07)) {
  # Call the fastglmPure function from fastglm package to fit the GLM
  fastglm::fastglmPure(
    x = x,
    y = y,
    weights = weights, # Observation weights
    family = family, # Family of the GLM
    offset = offset, # Offset for the linear predictor
    start = start, # Initial parameter estimates
    etastart = etastart, # Initial values for the linear predictor
    mustart = mustart, # Initial values for the mean
    method = 2L, # Method to be used
    tol = control$epsilon, # Tolerance for convergence
    maxit = control$maxit # Maximum number of iterations
  )
}

#' A Wrapper Function for arm::bayesglm.fit
#'
#' This function serves as a wrapper around arm::bayesglm.fit,
#' providing additional functionality and handling various input configurations.
#'
#' @importFrom arm bayesglm.fit
#'
#' @inheritParams arm::bayesglm.fit
#' @param control A list of control parameters including `maxit` for the maximum number of iterations (default is 100) and `epsilon` for the tolerance level (default is 1e-07).
#' @return A fitted generalized linear model with Bayesian prior distributions.
#'
.bayesglm_fit <- function(x,
                          y,
                          weights = rep(1, NROW(y)),
                          family = stats::gaussian(),
                          offset = rep(0, NROW(y)),
                          start = NULL,
                          mustart = NULL,
                          etastart = NULL,
                          control = list(maxit = 100, epsilon = 1e-07)) {
  arm::bayesglm.fit(
    x = x,
    y = y,
    weights = weights,
    start = start,
    etastart = etastart,
    mustart = mustart,
    offset = offset,
    family = family,
    control = control,
    prior.scale = 1,
    prior.mean.for.intercept = -14
  )
}

#' Fit a Generalized Linear Model (GLM) with Various Family Options
#'
#' This function provides a wrapper for various GLM fitting methods and returns coefficients
#' and other statistics of the fitted model.
#'
#' @param formula A formula object specifying the model. Default is "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))".
#' @param data A data frame containing the variables in the model.
#' @param family The family for the generalized linear model ("bayes.poisson", "poisson", "negative.binomial", "bayes.negative.binomial"). Defaults to "bayes.poisson".
#' @return A data table with the coefficients, p-values, and other statistics for each variable in the model.
#' @examples
#' # Your example code here, for example:
#' # your_formula <- "MutationNumber ~ isTarget + CNA + isTarget:CNA+ Mutation + offset(log(ntAtRisk))"
#' # fit_glm(formula = your_formula, data = your_data, family = "bayes.poisson")
#' @export
fit_glm <- function(formula = "MutationNumber ~ isTarget + CNA + isTarget:CNA + Mutation + offset(log(ntAtRisk))",
                    data,
                    family = c("bayes.poisson", "poisson", "negative.binomial", "bayes.negative.binomial")) {
  # Set the `family` to a validated value
  family <- base::match.arg(family)

  # Convert formula to formula object
  formula <- stats::as.formula(formula)

  # Try fitting the model using tryCatch to handle errors
  glm_model <- tryCatch(
    {
      if (family == "bayes.poisson") {
        arm::bayesglm(formula,
          data = data,
          family = stats::poisson(link = "log"),
          control = stats::glm.control(maxit = 100),
          prior.scale = 1
        ) %>% invisible()
      } else if (family == "poisson") {
        mm <- stats::model.matrix(formula, data = data)
        fastglm::fastglm(
          x = mm,
          y = data[["MutationNumber"]],
          family = stats::poisson(link = "log"),
          data = data,
          maxit = 100,
          method = 2L
        ) %>% invisible()
      } else if (family == "bayes.negative.binomial") {
        MASS::glm.nb(formula,
          data = data,
          method = ".bayesglm_fit"
        ) %>% invisible()
      } else if (family == "negative.binomial") {
        MASS::glm.nb(formula,
          data = data,
          control = stats::glm.control(maxit = 100),
          method = ".fastglm_fit"
        ) %>% invisible()
      }
    },
    error = function(err) {
      message("Error while fitting the model. If you are not using 'bayes.poisson' consider changing to it")
      NULL
    },
    finally = {}
  )

  # If model fitting fails, return a table with NAs
  if (is.null(glm_model)) {
    coefs <- data.table::data.table(
      coefName = NA_character_,
      coef = NA_real_,
      pVal = NA_real_,
      upperTailPval = NA_real_,
      lowerTailPval = NA_real_,
      Std.error = NA_real_,
      converged = FALSE
    )
    return(coefs)
  } else {
    if (!stringr::str_detect(family, "bayes")) {
      class(glm_model) <- "fastglm"
    }

    coefs <- summary(glm_model)$coefficients
    aCoef <- coefs[, 1L]
    aStd.error <- coefs[, 2L]
    aZ.Value <- coefs[, 3L]
    aPval <- coefs[, 4L]
    if (colnames(coefs)[3L] == "t value") {
      aupperTailPval <- stats::pt(coefs[, 3L], df = glm_model$df.residual, lower.tail = FALSE)
    } else {
      aupperTailPval <- stats::pnorm(coefs[, 3L], lower.tail = FALSE)
    }
    alowerTailPval <- 1 - aupperTailPval

    coefs <- data.table::data.table(
      coefName = rownames(coefs),
      coef = aCoef,
      pVal = aPval,
      upperTailPval = aupperTailPval,
      lowerTailPval = alowerTailPval,
      Std.error = aStd.error,
      converged = glm_model$converged
    )
    coefs[coefName == "(Intercept)", coefName := "Intercept"][]
    return(coefs)
  }
}
