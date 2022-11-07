#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics augment
#' @export
generics::augment

#' Tidy Model Results
#'
#' @param x A `fibre` model object
#' @param effects Which effects do you want tidied? One of: `"fixed"`, 
#' for fixed effects, `"random"` for random effects, or `"hyper"` for
#' the hyper-parameters of the random effects. Can also be `"rates"`,
#' which is a synonym for `"random"`, since the random effects are
#' rates of trait evolution along phylogenetic edges.
#' @param conf.type What kind of confidence interval. Choices are:
#' `"cred.int"` for approximate Bayesian marginal credible intervals. or
#' `"marginals"` for the full approximate marginal distributions, as a 
#' `data.frame` with `value` and `y.value` columns. `value` is the value
#' of the parameter, and `y.value` is the marginal posterior density 
#' (e.g. what `value` is the x axis and `y.value` is the y axis when plotting
#' the posterior density). 
#' @param indexes If `effects = "random"` or `effects = "rates"`, this is a 
#' vector of indices to retrieve particular random effects. Default is to return
#' all random effects, however, this can be slow for retrieving the marginals.
#' @param ... Not used.
#'
#' @return A tidy `tibble` with information about the fitted model parameters.
#' @export
tidy.fibre <- function(x, effects = c("fixed", "rates", "random", "hyper"),
                       conf.type = c("cred.int", "marginals"), 
                       indexes = NULL, ...) {
  
  effects <- match.arg(effects)
  conf.type <- match.arg(conf.type)
  
  if(conf.type == "cred.int") {
    tab <- switch(effects,
                  fixed = x$fixed %>%
                    dplyr::rename(conf.low = `0.025quant`,
                                  conf.high = `0.975quant`),
                  rates = ,
                  random = x$random %>%
                    purrr::imap_dfr(~ .x %>%
                                     dplyr::mutate(parameter = .y) %>%
                                      dplyr::rename(index = ID, conf.low = `0.025quant`,
                                                    conf.high = `0.975quant`) %>%
                                      dplyr::select(parameter, 
                                                    index,
                                                    mean,
                                                    sd,
                                                    conf.low,
                                                    conf.high)))
  } else {
    tab <- switch(effects,
                  fixed = tidy_marginal(x$model$marginals.fixed))
  }
  
  
}

tidy_marginal <- function(x) {
  purrr::imap_dfr(x,
                  ~ .x %>%
                    as.data.frame() %>%
                    dplyr::mutate(parameter = .y) %>%
                    dplyr::select(parameter, value = .data$x, y.value = .data$y))
}