#' Title
#'
#' @param formula formula specifying the response (LHS) and set of predictors (RHS) for the
#' phylogenetic model. The default \code{ ~ 1} is special shorthand specifying a model of
#' a response with only a phylogenetic model (e.g. no covariates), in which case, the \code{data}
#' argument must contain a vector or column matrix/data.frame with the response variable(s)
#' (e.g. species traits).
#' @param phy A \code{phylo} object containing the phylogeny to be used for the phylogenetic trait
#' model.
#' @param data A matrix or data.frame (or vector) containing the variables referred to in
#' \code{formula} argument.
#' @param phy_match This argument specifies how the data should be matched to the phylogeny. The
#' default \code{"auto"}, will match by names or rownames if they exist in \code{data}, or match
#' by order if no names are present. You can also explicitly specify how to match (recommended).
#' Option include "names", to match by names (throwing an error if they aren't present), "order", to
#' match by ordering (element or rows in \code{data} are in the same order as \code{phy$tip.label}),
#' or any length 1 character vector that refers to a column name in \code{data}, containing a character
#' vector. In this case, this character vector will be used to match the data to the tip labels in
#' \code{phy}. This is the recommended best way, as it is the most compatible with "tidy" principles,
#' making sure species names are treated as data and stay associated with related data.
#' @param family Family for the distribution of errors. Default is "gaussian", but can be any
#' supported by \code{inla}. See \code{names(inla.models()$likelihood)} for a list of possibilities
#' and use \code{\link[INLA]{inla.doc}} for details of all available families.
#' @param rate_model The model to use for the evolutionary rate variation along the phylogeny
#' (in other words, the choice of prior for rates). Current choices are:
#' \itemize{
#' \item{bayes_ridge}{Rates are completely independent along all branched but are shrunk
#' towards zero by an independent Normal prior. This is classic Bayesian Ridge Regression and has
#' a single hyperparameter whose prior determines the degree of shrinkage.}
#' \item{temporal_rates}{Rates are constrained to be similar if they are close together
#' in time. Has a single hyperparameter whose prior determines the degree of temporal
#' "smoothing". Heavy smoothing forces rates to change slowly through time, light smoothing
#' allows rates to change quickly (e.g. they can be "wiggly"). Specify a particular model of
#' temporal autocorrelation with the \code{temporal_model} argument.}
#' }
#' @param temporal_model
#' @param fit
#' @param aces
#' @param hyper_priors
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
phinla <- function(formula = ~ 1, phy, data = NULL,
                   phy_match = "auto",
                   family = "gaussian",
                   rate_model = c("bayes_ridge", "temporal_rates",
                             "brownian_rates", "temporal_plus_brownian"),
                   temporal_model = "crw1",
                   fit = TRUE, aces = TRUE,
                   hyper_priors = NULL,
                   ...) {

  rate_model <- match.arg(rate_model)

  phy_mat <- RRphylo::makeL(phy)[ , -1]

  node_names <- colnames(phy_mat)
  tip_names <- phy$tip.label

  if(phy_match == "auto") {
    if(is.null(dim(data))) {
      if(is.null(names(data))) {
        data <- dplyr::tibble(node_name = tip_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = phy$tip.label,
                                         y = data))

      } else {
        data <- dplyr::tibble(node_name = tip_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = names(data),
                                         y = data))
      }

    } else {
      if(is.null(rownames(data))) {
        data <- dplyr::tibble(node_name = tip_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = phy$tip.label) %>%
                             dplyr::bind_cols(as.data.frame(data)))
      } else {
        data <- dplyr::tibble(node_name = tip_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = rownames(data)) %>%
                             dplyr::bind_cols(as.data.frame(data)))
      }
    }
  } else {
    if(phy_match == "names") {
      if(is.null(dim(data))) {
        data <- dplyr::tibble(node_name = tip_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = names(data),
                                         y = data))
      } else {
        data <- dplyr::tibble(node_name = tip_names) %>%
          dplyr::left_join(dplyr::tibble(node_name = rownames(data)) %>%
                             dplyr::bind_cols(as.data.frame(data)))
      }
    } else {
      if(phy_match == "order") {
        if(is.null(dim(data))) {
          data <- dplyr::tibble(node_name = tip_names) %>%
            dplyr::left_join(dplyr::tibble(node_name = phy$tip.label,
                                           y = data))
        } else {
          data <- dplyr::tibble(node_name = tip_names) %>%
            dplyr::left_join(dplyr::tibble(node_name = phy$tip.label) %>%
                               dplyr::bind_cols(as.data.frame(data)))
        }
      } else {
        if(!is.null(dim(data))) {
          data <- as.data.frame(data)
          if(!phy_match %in% colnames(data)) {
            stop("Name provided in phy_match does not match any column names in data.")
          } else {
            data <- data %>%
              dplyr::rename(node_name = tip_names)
            data <- dplyr::tibble(node_name = tip_names) %>%
              dplyr::left_join(data)
          }
        } else {
          stop('If phy_match does not equal "auto", "names", or "order", then data must
               be a matrix or data.frame, not a vector')
        }
      }
    }
  }

  te <- terms(formula)
  if(attr(te, "response") == 0) {
    formula <- update(formula, y ~ .)
  }

  dat <- model.frame(formula, data, na.action = "na.pass") %>%
    as.data.frame()
  nam <- data$node_name

  dat <- dat %>%
    dplyr::mutate(`root` = 1)


  tip_indexes <- 1:nrow(dat)
  node_indexes <- 1:ncol(phy_mat)

  A_mat <- Matrix::Matrix(phy_mat[nam, ])

  resp <- all.vars(formula[[2]])

  phy_stack <- INLA::inla.stack(data = list(y = dat[ , resp]),
                              A = list(A_mat, 1),
                              effects = list(node_id = node_indexes,
                                             root = dat$root),
                              tag = "rates")

  if(aces) {
    aces_A_mat <- RRphylo::makeL1(phy)[ , -1]
    tip_mat <- matrix(0, nrow = nrow(aces_A_mat), ncol = length(phy$tip.label))
    colnames(tip_mat) <- phy$tip.label
    aces_A_mat <- cbind(aces_A_mat, tip_mat)
    aces_A_mat <- aces_A_mat[ , colnames(phy_mat)]
    aces_A_mat <- Matrix::Matrix(aces_A_mat)

    aces_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(aces_A_mat))),
                                   A = list(aces_A_mat, 1),
                                   effects = list(node_id = 1:ncol(aces_A_mat),
                                                  root = rep(1, nrow(aces_A_mat))),
                                   tag = "aces")

    full_stack <- INLA::inla.stack(phy_stack, aces_stack)
  }

  pc.prec = list(prec = list(prior = "pc.prec", param = c(5, 0.01)))

  if(rate_model == "bayes_ridge") {
    inla_form <- y ~ 0 + root + f(node_id, model = "iid",
                                  constr = FALSE,
                                  hyper = pc.prec)

  }

  if(aces) {
    fit_modes <- INLA::inla(inla_form,
                       data = INLA::inla.stack.data(phy_stack),
                       family = family,
                       control.predictor = list(A = inla.stack.A(phy_stack),
                                                compute = FALSE))

    fit <- INLA::inla(inla_form,
                       data = INLA::inla.stack.data(full_stack),
                       family = family,
                       control.predictor = list(A = inla.stack.A(full_stack),
                                                compute = TRUE),
                       control.mode = list(theta = fit_modes$mode$theta, restart = FALSE))
  } else {
    fit <- INLA::inla(inla_form,
                      data = INLA::inla.stack.data(phy_stack),
                      family = family,
                      control.predictor = list(A = inla.stack.A(phy_stack),
                                               compute = TRUE))
  }

  fit

}
