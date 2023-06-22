#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @export
autoplot.fibre <- function(object, element = c("params", "hyper", "fixed", "rates", "aces"), param_names = NULL, n = 10000, suppress_intercept = FALSE, scale = NULL, standardise_height = TRUE, ...) {
  element <- match.arg(element)
  switch(element,
         params = plot_params(object, n = n, suppress_intercept = suppress_intercept, scale = scale, param_names = param_names, standardise_height = standardise_height))
}

plot_params <- function(object, n = 10000, suppress_intercept = FALSE, scale = NULL, param_names = param_names, standardise_height = TRUE) {

  rlang::check_installed(c("ggridges", "INLA", "ggforce"))

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

  if(!is.null(param_names)) {
    if(length(param_names) != 2) {
      rlang::abort("param_names should be length two, the first element referring to the random effects, the second to the fixed effects")
    }

    random_df <- dplyr::tibble(parameter = param_names[[1]]) %>%
      left_join(object$hyper)

    if(length(param_names[[2]]) > 1) {
      split_fixed <- TRUE
      fixed_df <- purrr::map2(param_names[[2]], names(param_names[[2]]),
                              ~ dplyr::tibble(parameter = .x) %>%
                                dplyr::left_join(object$fixed) %>%
                                dplyr::mutate(effect_type = paste("Fixed Effect:", .y)))
      # effect_types <- fixed_df$effect_type
      # fixed_df <- fixed_df %>%
      #   dplyr::select(-any_of("effect_type"))
    } else {
      fixed_df <- dplyr::tibble(parameter = param_names[[2]]) %>%
        left_join(object$fixed)
      split_fixed <- FALSE
    }

    if(!is.null(names(param_names[[1]]))) {
      random_df$parameter <- names(param_names[[1]])
    }
    if(split_fixed) {

      for(i in seq_along(param_names[[2]])) {
        if(!is.null(names(param_names[[2]][[i]]))) {
          fixed_df[[i]]$parameter <- names(param_names[[2]][[i]])
        }
      }

    } else {

      if(!is.null(names(param_names[[2]]))) {
        fixed_df$parameter <- names(param_names[[2]])
      }

    }
  } else {

    random_df <- object$hyper
    fixed_df <- object$fixed
    split_fixed <- FALSE

  }

  if(split_fixed) {
    for(i in seq_along(fixed_df)) {
      fixed_df[[i]]$parameter <- factor(fixed_df[[i]]$parameter, levels = fixed_df[[i]]$parameter)
    }

  } else {

    fixed_df$parameter <- factor(fixed_df$parameter, levels = fixed_df$parameter)

  }

  random_df$parameter <- factor(random_df$parameter, levels = random_df$parameter)

  if(split_fixed) {
    fixed_samps <- purrr::map_dfr(fixed_df,
                              ~ get_marg_samps(.x$marginal, .x$parameter, n = n) %>%
                                dplyr::mutate(effect_type = .x$effect_type[1]))
  } else {
    fixed_samps <- get_marg_samps(fixed_df$marginal, fixed_df$parameter, n = n) %>%
      dplyr::mutate(effect_type = "Fixed Effects")
  }

  hyperpar_samps <- get_marg_samps(random_df$marginal, random_df$parameter, n = n) %>%
    dplyr::mutate(effect_type = "Random Effects")

  if(split_fixed) {
    fixed_df <- dplyr::bind_rows(fixed_df)
    fixed_effect_types <- unique(fixed_df$effect_type)
  } else {
    fixed_effect_types <- "Fixed Effects"
  }

  all_samps <- dplyr::bind_rows(fixed_samps,
                                hyperpar_samps) %>%
    dplyr::mutate(effect_type = factor(effect_type,
                                       levels = c("Random Effects", fixed_effect_types)))

  if(split_fixed) {

    ci <- bind_rows(fixed_df %>%
                    dplyr::select(var = parameter,
                                  lower = `0.025quant`,
                                  upper = `0.975quant`,
                                  mean,
                                  effect_type),
                    random_df %>%
                      dplyr::select(var = parameter,
                                    lower = `0.025quant`,
                                    upper = `0.975quant`,
                                    mean) %>%
                      dplyr::mutate(effect_type = "Random Effects")) %>%
      dplyr::mutate(effect_type = factor(effect_type,
                                         levels = levels(all_samps$effect_type)))

  } else {
    ci <- bind_rows(fixed_df %>%
                    dplyr::select(var = parameter,
                                  lower = `0.025quant`,
                                  upper = `0.975quant`,
                                  mean) %>%
                      dplyr::mutate(effect_type = "Fixed Effects"),
                    random_df %>%
                      dplyr::select(var = parameter,
                                    lower = `0.025quant`,
                                    upper = `0.975quant`,
                                    mean) %>%
                      dplyr::mutate(effect_type = "Random Effects")) %>%
      dplyr::mutate(effect_type = factor(effect_type,
                                         levels = levels(all_samps$effect_type)))
  }



  sig_vars <- ci %>%
    dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
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
                                       levels = levels(all_samps$effect_type)))

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
                                         levels = levels(all_samps$effect_type)))
  }

  pal <- c("#fc8d62", "#8da0cb")
  if(standardise_height) {
    p <- ggplot2::ggplot(samps, ggplot2::aes(val, var, group = var, height = after_stat(ndensity)))
  } else {
    p <- ggplot2::ggplot(samps, ggplot2::aes(val, var, group = var, height = after_stat(density)))
  }

  if(scale_multi) {
    # p <- p +
    #   ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig, scale = scale),
    #                                 stat = "density", adjust = 2, color = "gray70",
    #                                 rel_min_height = 0.01)

    p <- p +
      ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig, scale = scale),
                                    stat = "density", adjust = 2, color = "gray70",
                                    rel_min_height = 0.01)
  } else {
    # p <- p +
    #   ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig),
    #                                 stat = "density", adjust = 2, color = "gray70",
    #                                 scale = scale,
    #                                 rel_min_height = 0.01)

    p <- p +
      ggridges::geom_density_ridges(ggplot2::aes(alpha = sig, fill = sig),
                                    stat = "density", adjust = 2, color = "gray70",
                                    scale = scale,
                                    rel_min_height = 0.01)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(x = mean, y = var), data = ci, inherit.aes = FALSE) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper, y = var), data = ci,
                   inherit.aes = FALSE, height = 0.1) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    scale_y_discrete(expand = c(0.01, 0)) +
    ggplot2::scale_alpha_manual(values = c(0.8, 0.2)) +
    ggplot2::scale_fill_manual(values = rev(pal)) +
    ggforce::facet_col(ggplot2::vars(effect_type),
                       scales = "free", space = "free") +
    ggplot2::ylab("") +
    ggplot2::xlab("Estimate") +
    ggridges::theme_ridges() +
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
                        values_to = "val") %>%
    dplyr::group_by(var) %>%
    filter(val < quantile(val, 0.999) & val > quantile(val, 0.001)) %>%
    ungroup()

}


