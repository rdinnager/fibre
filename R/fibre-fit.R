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
fibre.data.frame <- function(x, y, 
                             intercept = TRUE, 
                             engine = c("inla", "glmnet", "mgcv"),
                             engine_options = list(),
                             ...) {
  engine <- match.arg(engine)
  processed <- hardhat::mold(x, y)
  fibre_bridge(processed, engine, engine_options, ...)
}

# XY method - matrix

#' @export
#' @rdname fibre
fibre.matrix <- function(x, y, 
                         intercept = TRUE,
                         engine = c("inla", "glmnet", "mgcv"),
                         engine_options = list(),
                         ...) {
  engine <- match.arg(engine)
  processed <- hardhat::mold(x, y)
  fibre_bridge(processed, engine, engine_options, ...)
}

# Formula method

#' @export
#' @rdname fibre
fibre.formula <- function(formula, data, 
                          intercept = TRUE,
                          family = "gaussian",
                          engine = c("inla", "glmnet", "mgcv"),
                          engine_options = list(),
                          ...) {
  engine <- match.arg(engine)
  processed <- hardhat::mold(formula, data, 
                             blueprint = fibre_formula_blueprint(intercept = intercept))
  
  if(engine == "glmnet" && length(processed$extras$model_info) > 1) {
    rlang::abort('engine = "glmnet" currently only supports a single bre() call in a model.')
  }
  fibre_bridge(processed, family, engine, engine_options, ...)
}

# Recipe method

#' @export
#' @rdname fibre
fibre.recipe <- function(x, data, 
                         intercept = TRUE, 
                         engine = c("inla", "glmnet", "mgcv"),
                         engine_options = list(),
                         ...) {
  engine <- match.arg(engine)
  processed <- hardhat::mold(x, data)
  fibre_bridge(processed, engine, engine_options, ...)
}

# ------------------------------------------------------------------------------
# Bridge

fibre_bridge <- function(processed, family, engine, engine_options, ...) {
  
  predictors <- processed$predictors
  outcomes <- processed$outcomes
  offset <- processed$extras$offset
  pfcs <- purrr::map(processed$extras$model_info,
                     "phyf")
  rate_dists <- purrr::map(processed$extras$model_info,
                           "rate_dist")
  hypers <- purrr::map(processed$extras$model_info,
                       "hyper")
  latents <- purrr::map(processed$extras$model_info,
                        "latent")
  mixture_ofs <- purrr::map(processed$extras$model_info,
                           "mixture_of")

  
  fit <- fibre_impl(predictors, outcomes,
                    offset, pfcs,
                    rate_dists, hypers,
                    latents, family,
                    engine,
                    engine_options)
  
  #return(fit)
  
  switch(engine,
         inla = fibre_process_fit_inla(fit, processed$blueprint,
                                       predictors,
                                       pfcs,
                                       rate_dists))

}


# ------------------------------------------------------------------------------
# Implementation

fibre_impl <- function(predictors, outcomes,
                       offset, pfcs,
                       rate_dists, hypers,
                       latents,
                       family,
                       engine,
                       engine_options) {
  
  if(engine == "glmnet") {
    rlang::abort('engine = "glmnet" does not support argument latent > 0')
  }
  
  dat_list <- switch(engine,
                     inla = shape_data_inla(pfcs,
                                            predictors,
                                            outcomes,
                                            latents),
                     glmnet = shape_data_glmnet(pfcs,
                                                predictors,
                                                outcomes))
  
  form <- switch(engine, 
                 inla = make_inla_formula(dat_list$dat, dat_list$y),
                 glmnet = NULL)
  
  
  family <- get_families(family, colnames(dat_list$y))

  
  if(engine == "inla") {
    #names(dat_list$y) <- backtick_names(names(dat_list$y))
    #names(dat_list$dat) <- backtick_names(names(dat_list$dat))
    hypers <- form$hypers
    form <- form$form
    inla_dat <- INLA::inla.stack(data = list(y = dat_list$y),
                                 A = list(dat_list$A),
                                 effects = list(dat_list$dat),
                                 compress = FALSE,
                                 remove.unused = FALSE)
    inla_options <- list(control.predictor = list(A = INLA::inla.stack.A(inla_dat),
                                                  link = 1,
                                                  compute = TRUE),
                         inla.mode = "experimental")
    inla_options <- utils::modifyList(inla_options, engine_options, keep.null = TRUE)
    
    n_re <- length(hypers$re) 
    if(n_re > 0) {
      hypers_re <- hypers[seq_along(hypers$re)]
      names(hypers_re) <- hypers$re
      rlang::env_bind(rlang::f_env(form), !!!hypers_re)
    } 
    
    if(length(hypers$latent) > 0) {
      hypers_latent <- hypers[seq_along(hypers$latent) + n_re]
      names(hypers_latent) <- hypers$latent
      hypers_copy <- rep(list(list(beta = list(fixed = FALSE))), length(hypers$copy))
      names(hypers_copy) <- hypers$copy
      rlang::env_bind(rlang::f_env(form), !!!hypers_latent)
      rlang::env_bind(rlang::f_env(form), !!!hypers_copy)
    }
    
  }
  
  #return(list(inla_dat, form))
  
  fit <- switch(engine,
                inla = rlang::exec(INLA::inla, formula = form,
                                   data = INLA::inla.stack.data(inla_dat),
                                   family = family,
                                   !!!inla_options),
                rlang::abort("Invalid engine argument")
  )
  
  list(fit = fit, renamer = dat_list$renamer)
  
}
