phinla <- function(formula, phy, data = NULL,
                   family = "gaussian",
                   rate_model = c("bayes_ridge", "temporal_rates",
                             "brownian_rates"),
                   fit = TRUE, id_col = NULL) {

  rate_model <- match.arg(rate_model)

  phy_mat <- RRphylo::makeL(phy)
  node_names <- colnames(phy_mat)

  if(is.null(id_col)) {
    if(is.null(dim(data))) {
      if(is.null(names(data))) {
        data <- dplyr::tibble(node_name = node_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = phy$tip.label,
                                         y = data))

      } else {
        data <- dplyr::tibble(node_name = node_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = names(data),
                                         y = data))
      }

    } else {
      if(is.null(rownames(data))) {
        data <- dplyr::tibble(node_name = node_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = phy$tip.label) %>%
                             dplyr::bind_cols(as.data.frame(data)))
      } else {
        dplyr::tibble(node_name = node_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = rownames(data)) %>%
                             dplyr::bind_cols(as.data.frame(data)))
      }
    }
  } else {
    if(!inherits(data, "data.frame")) {
      stop("You have specified an id column but data is not a data.frame.")
    } else {
      colnames(data)[colnames(data) == id_col] <- "tip_name"
    }
  }

  dat <- model.frame(formula, data)

}
