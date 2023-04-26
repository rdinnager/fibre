#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @export
autoplot.fibre <- function(object, element = c("params", "hyper", "fixed", "rates", "aces"), n = 10000, suppress_intercept = FALSE, scale = NULL, ...) {
  element <- match.arg(element)
  switch(element,
         params = plot_params(object, n = n, suppress_intercept = suppress_intercept, scale = scale))
}

plot_params <- function(object, n = 10000, suppress_intercept = FALSE, scale = NULL) {

  rlang::check_installed(c("ggridges", "INLA"))

  if(!is.null(scale)) {
    if(length(scale) > 1) {
      scale_multi <- TRUE
      if(is.null(names(scale))) {
        names(scale) <- c("Random Effects", "Fixed Effects")
      }
      scale <- dplyr::tibble(effect_type = names(scale),
                             scale = scale)
    } else {
      scale_multi <- FALSE
    }
  } else {
    scale <- 1
  }

  fixed_samps <- get_marg_samps(object$fixed$marginal, object$fixed$parameter, n = n) %>%
    mutate(effect_type = "Fixed Effects")
  hyperpar_samps <- get_marg_samps(object$hyper$marginal, object$hyper$parameter, n = n) %>%
    mutate(effect_type = "Random Effects")

  all_samps <- dplyr::bind_rows(fixed_samps,
                                hyperpar_samps) %>%
    dplyr::mutate(effect_type = factor(effect_type,
                                       levels = c("Random Effects", "Fixed Effects")))

  ci <- bind_rows(object$fixed %>%
                    dplyr::select(var = parameter,
                                  lower = `0.025quant`,
                                  upper = `0.975quant`,
                                  mean) %>%
                    dplyr::mutate(effect_type = "Fixed Effects"),
                  object$hyper %>%
                    dplyr::select(var = parameter,
                                  lower = `0.025quant`,
                                  upper = `0.975quant`,
                                  mean) %>%
                    dplyr::mutate(effect_type = "Random Effects")) %>%
    dplyr::mutate(effect_type = factor(effect_type,
                                       levels = c("Random Effects", "Fixed Effects")))
  # ci <- all_samps %>%
  #   dplyr::group_by(var, effect_type) %>%
  #   dplyr::summarise(lower = quantile(val, 0.025),
  #                    upper = quantile(val, 0.975),
  #                    mean = mean(val),
  #                    .groups = "drop_last")

  sig_vars <- ci %>%
    dplyr::mutate(sig = ifelse(effect_type == "random",
                               "sig",
                               ifelse(sign(lower) == sign(upper),
                                      "sig",
                                      "no_sig"))) %>%
    dplyr::select(var, sig)

  samps <- all_samps %>%
    dplyr::left_join(sig_vars, by = "var") %>%
    dplyr::group_by(var) %>%
    dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(effect_type = factor(effect_type,
                                       levels = c("Random Effects", "Fixed Effects")))

  if(suppress_intercept) {

    samps <- samps %>%
      filter(!grepl("(Intercept)", var))

    ci <- ci %>%
      filter(!grepl("(Intercept)", var))

  }

  if(scale_multi) {
    samps <- samps %>%
      left_join(scale) %>%
      dplyr::mutate(effect_type = factor(effect_type,
                                         levels = c("Random Effects", "Fixed Effects")))
  }

  pal <- c("#fc8d62", "#8da0cb")
  p <- ggplot2::ggplot(samps, ggplot2::aes(val, var, height = after_stat(density)))

  if(scale_multi) {
    p <- p +
      ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig, scale = scale),
                                    stat = "density", adjust = 2, color = "gray70")
  } else {
    p <- p +
      ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig),
                                    stat = "density", adjust = 2, color = "gray70",
                                    scale = scale)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(x = mean, y = var), data = ci, inherit.aes = FALSE) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper, y = var), data = ci,
                   inherit.aes = FALSE, height = 0.1) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::scale_alpha_manual(values = c(0.8, 0.2)) +
    ggplot2::scale_fill_manual(values = rev(pal)) +
    ggplot2::facet_wrap(ggplot2::vars(effect_type), nrow = 2, scales = "free") +
    ggplot2::ylab("") +
    ggplot2::xlab("Estimate") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
          axis.text = ggplot2::element_text(size = 14),
          strip.text = ggplot2::element_text(size = 16))

  p

}

get_marg_samps <- function(marg, names, n = 10000) {

  purrr::map(marg,
             ~ INLA::inla.rmarginal(n, .x)) %>%
    setNames(names) %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "var",
                        values_to = "val")

}


