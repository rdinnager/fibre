#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics augment
#' @export
generics::augment

tidy.fibre <- function(x, effects = c("fixed", "rates", "random", "hyper"),
                       conf.type = c("cred.int", "marginals"), ...) {
  
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