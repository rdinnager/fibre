#' Title
#'
#' @param x A fitted model object produced by \code{\link{phybrr}}
#' @param type What kind of posterior summary to return?
#' @param n If \code{type = "samples"}, how many samples to return?
#' @param p If \code{type = "hpd"}, what alpha levels to use?
#'
#' @return For all types except "hpd", "ci", and "marginals", a
#' numeric vector, otherwise a list for "hpd" and "marginals", and a matrix
#' for "ci".
#' @export
get_rates <- function(x, type = c("marginals", "samples",
                                  "mode", "mean", "median",
                                  "lower", "upper", "ci", "hpd",
                                  "sd"),
                      n = 1,
                      p = 0.05) {

  type <- match.arg(type)

  rate_ind <- attr(x, "indexes")$rates
  switch(type,
         marginals = x$marginals.fitted.values[rate_ind],
         samples = lapply(rate_ind,
                          function(y) INLA::inla.rmarginal(x$marginals.fitted.values[[y]], n = n)),
         mode = x$summary.fitted.values$mode[rate_ind],
         mean = x$summary.fitted.values$mean[rate_ind],
         median = x$summary.fitted.values$`0.5quant`[rate_ind],
         lower = x$summary.fitted.values$`0.025quant`[rate_ind],
         upper = x$summary.fitted.values$`0.975quant`[rate_ind],
         ci = cbind(x$summary.fitted.values$`0.025quant`[rate_ind],
                x$summary.fitted.values$`0.975quant`[rate_ind]),
         hpd = lapply(rate_ind,
                      function(y) INLA::inla.hpdmarginal(p, x$marginals.fitted.values[[y]])),
         sd = x$summary.fitted.values$sd[rate_ind]
         )
}

#' Title
#'
#' @param x A fitted model object produced by \code{\link{phybrr}}
#' @param type What kind of posterior summary to return?
#' @param n If \code{type = "samples"}, how many samples to return?
#' @param p If \code{type = "hpd"}, what alpha levels to use?
#'
#' @return For all types except "hpd", "ci", and "marginals", a
#' numeric vector, otherwise a list for "hpd" and "marginals", and a matrix
#' for "ci".
#' @export
get_aces <- function(x, type = c("marginals", "samples",
                                 "mode", "mean", "median",
                                 "lower", "upper", "ci", "hpd",
                                 "sd"),
                      n = 1,
                      p = 0.05) {

  type <- match.arg(type)

  rate_ind <- attr(x, "indexes")$aces
  switch(type,
         marginals = x$marginals.fitted.values[rate_ind],
         samples = lapply(rate_ind,
                          function(y) INLA::inla.rmarginal(x$marginals.fitted.values[[y]], n = n)),
         mode = x$summary.fitted.values$mode[rate_ind],
         mean = x$summary.fitted.values$mean[rate_ind],
         median = x$summary.fitted.values$`0.5quant`[rate_ind],
         lower = x$summary.fitted.values$`0.025quant`[rate_ind],
         upper = x$summary.fitted.values$`0.975quant`[rate_ind],
         ci = cbind(x$summary.fitted.values$`0.025quant`[rate_ind],
                    x$summary.fitted.values$`0.975quant`[rate_ind]),
         hpd = lapply(rate_ind,
                      function(y) INLA::inla.hpdmarginal(p, x$marginals.fitted.values[[y]])),
         sd = x$summary.fitted.values$sd[rate_ind]
  )
}

#' Title
#'
#' @param x A fitted model object produced by \code{\link{phybrr}}
#' @param type What kind of posterior summary to return?
#' @param n If \code{type = "samples"}, how many samples to return?
#' @param p If \code{type = "hpd"}, what alpha levels to use?
#'
#' @return For all types except "hpd", "ci", and "marginals", a
#' numeric vector, otherwise a list for "hpd" and "marginals", and a matrix
#' for "ci".
#' @export
get_tces <- function(x, type = c("marginals", "samples",
                                 "mode", "mean", "median",
                                 "lower", "upper", "ci", "hpd",
                                 "sd"),
                      n = 1,
                      p = 0.05) {

  type <- match.arg(type)

  rate_ind <- attr(x, "indexes")$tces
  switch(type,
         marginals = x$marginals.fitted.values[rate_ind],
         samples = lapply(rate_ind,
                          function(y) INLA::inla.rmarginal(x$marginals.fitted.values[[y]], n = n)),
         mode = x$summary.fitted.values$mode[rate_ind],
         mean = x$summary.fitted.values$mean[rate_ind],
         median = x$summary.fitted.values$`0.5quant`[rate_ind],
         lower = x$summary.fitted.values$`0.025quant`[rate_ind],
         upper = x$summary.fitted.values$`0.975quant`[rate_ind],
         ci = cbind(x$summary.fitted.values$`0.025quant`[rate_ind],
                    x$summary.fitted.values$`0.975quant`[rate_ind]),
         hpd = lapply(rate_ind,
                      function(y) INLA::inla.hpdmarginal(p, x$marginals.fitted.values[[y]])),
         sd = x$summary.fitted.values$sd[rate_ind]
  )
}


#' Title
#'
#' @param x A fitted model object produced by \code{\link{phybrr}}
#' @param type What kind of posterior summary to return?
#' @param n If \code{type = "samples"}, how many samples to return?
#' @param p If \code{type = "hpd"}, what alpha levels to use?
#'
#' @return For all types except "hpd", "ci", and "marginals", a
#' numeric vector, otherwise a list for "hpd" and "marginals", and a matrix
#' for "ci".
#' @export
get_tips <- function(x, type = c("marginals", "samples",
                                 "mode", "mean", "median",
                                 "lower", "upper", "ci", "hpd",
                                 "sd"),
                     n = 1,
                     p = 0.05) {

  type <- match.arg(type)

  rate_ind <- attr(x, "indexes")$tips
  switch(type,
         marginals = x$marginals.fitted.values[rate_ind],
         samples = lapply(rate_ind,
                          function(y) INLA::inla.rmarginal(x$marginals.fitted.values[[y]], n = n)),
         mode = x$summary.fitted.values$mode[rate_ind],
         mean = x$summary.fitted.values$mean[rate_ind],
         median = x$summary.fitted.values$`0.5quant`[rate_ind],
         lower = x$summary.fitted.values$`0.025quant`[rate_ind],
         upper = x$summary.fitted.values$`0.975quant`[rate_ind],
         ci = cbind(x$summary.fitted.values$`0.025quant`[rate_ind],
                    x$summary.fitted.values$`0.975quant`[rate_ind]),
         hpd = lapply(rate_ind,
                      function(y) INLA::inla.hpdmarginal(p, x$marginals.fitted.values[[y]])),
         sd = x$summary.fitted.values$sd[rate_ind]
  )
}
