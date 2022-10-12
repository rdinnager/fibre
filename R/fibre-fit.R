#' Fit a `fibre`
#'
#' `fibre()` fits a model.
#'
#' @param x Depending on the context:
#'
#'   * A __data frame__ of predictors.
#'   * A __matrix__ of predictors.
#'   * A __recipe__ specifying a set of preprocessing steps
#'     created from [recipes::recipe()].
#'
#' @param y When `x` is a __data frame__ or __matrix__, `y` is the outcome
#' specified as:
#'
#'   * A __data frame__ with 1 numeric column.
#'   * A __matrix__ with 1 numeric column.
#'   * A numeric __vector__.
#'
#' @param data When a __recipe__ or __formula__ is used, `data` is specified as:
#'
#'   * A __data frame__ containing both the predictors and the outcome.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side.
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @return
#'
#' A `fibre` object.
#'
#' @examples
#' predictors <- mtcars[, -1]
#' outcome <- mtcars[, 1]
#'
#' # XY interface
#' mod <- fibre(predictors, outcome)
#'
#' # Formula interface
#' mod2 <- fibre(mpg ~ ., mtcars)
#'
#' # Recipes interface
#' library(recipes)
#' rec <- recipe(mpg ~ ., mtcars)
#' rec <- step_log(rec, disp)
#' mod3 <- fibre(rec, mtcars)
#'
#' @export
fibre <- function(x, ...) {
  UseMethod("fibre")
}

#' @export
#' @rdname fibre
fibre.default <- function(x, ...) {
  stop("`fibre()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# XY method - data frame

#' @export
#' @rdname fibre
fibre.data.frame <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  fibre_bridge(processed, ...)
}

# XY method - matrix

#' @export
#' @rdname fibre
fibre.matrix <- function(x, y, ...) {
  processed <- hardhat::mold(x, y)
  fibre_bridge(processed, ...)
}

# Formula method

#' @export
#' @rdname fibre
fibre.formula <- function(formula, data, ...) {
  processed <- hardhat::mold(formula, data)
  fibre_bridge(processed, ...)
}

# Recipe method

#' @export
#' @rdname fibre
fibre.recipe <- function(x, data, ...) {
  processed <- hardhat::mold(x, data)
  fibre_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

fibre_bridge <- function(processed, ...) {
  predictors <- processed$predictors
  outcome <- processed$outcomes[[1]]

  fit <- fibre_impl(predictors, outcome)

  new_fibre(
    coefs = fit$coefs,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation

fibre_impl <- function(predictors, outcome) {
  list(coefs = 1)
}