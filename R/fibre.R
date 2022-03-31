#' Phylogenetic Branch Regression (fibre)
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
#' \item{`"ridge"`}{Rates are completely independent along all branches but are shrunk
#' towards zero by an independent Normal prior. This is classic Bayesian Ridge Regression and has
#' a single hyperparameter whose prior determines the degree of shrinkage.}
#' \item{`"temporal"`}{Rates are constrained to be similar if they are close together
#' in 'time' (where time is defined as the distance from the root on the phylogeny).
#' Has a single hyperparameter whose prior determines the degree of temporal
#' "smoothing". Heavy smoothing forces rates to change slowly through time, light smoothing
#' allows rates to change quickly (e.g. they can be "wiggly"). Specify a particular model of
#' temporal autocorrelation with the \code{temporal_model} argument.}
#' }
#' @param fit
#' @param aces
#' @param hyper_priors
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fibre <- function(formula = ~ 1, phy, data = NULL,
                   phy_match = "auto",
                   family = "gaussian",
                   rate_model = "ridge",
                   rate_order = c("first_order", "second_order"),
                   fit = TRUE, aces = TRUE,
                   hyper = NULL,
                   obs_error = c("est", "one", "zero"),
                   ...) {

  rate_model <- match.arg(rate_model,
                          c("ridge", "temporal", "phylogenetic"),
                          several.ok = TRUE)
  obs_error <- match.arg(obs_error)
  rate_order <- match.arg(rate_order)

  fam_cont <- switch(obs_error,
                     est = list(),
                     one = list(hyper = list(prec = list(prior = "gaussian",
                                                         initial = 1,
                                                         fixed = TRUE))),
                     zero = list(hyper = list(prec = list(prior = "gaussian",
                                                          initial = 0,
                                                          fixed = TRUE))))

  message("Assembling model data and structure...")


  if(inherits(phy, "phylo")) {
    message("Generating root-to-tip matrix...")
    if(aces) {
      phy_mat <- make_root2tip(phy, return_nodes = "both",
                               return_type = "list",
                               sparse = TRUE,
                               order = rate_order,
                               return_ages = "temporal" %in% rate_model)
      if("temporal" %in% rate_model) {
        ages <- attr(phy_mat, "ages")
      }
      aces_A_mat <- phy_mat[[2]]
      phy_mat <- phy_mat[[1]]
    } else {
      phy_mat <- make_root2tip(phy, return_nodes = "tips",
                               sparse = TRUE,
                               order = rate_order,
                               return_ages = "temporal" %in% rate_model)
      if("temporal" %in% rate_model) {
        ages <- attr(phy_mat, "ages")
      }
    }
  } else {
    phy_mat <- phy[[2]]
    phy <- phy[[1]]
    if(aces) {
      aces_A_mat <- phy[[3]]
    }
  }

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
              dplyr::rename(node_name = dplyr::all_of(phy_match))
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
    vars <- setdiff(colnames(data), "node_name")
    if(length(vars) > 1) {
      formula <- update(formula, as.formula(paste0("cbind(", paste(vars, collapse = ", "),
                                   ") ~ .")))
    } else {
      formula <- update(formula, y ~ .)
    }
  }

  dat <- model.frame(formula, data, na.action = "na.pass")

  if(is.matrix(dat[ , 1])) {
    fits <- pbapply::pblapply(as.data.frame(dat[ , 1]), function(k) {
      names(k) <- data$node_name
      suppressMessages(fibrer(formula = ~ 1, phy = c(list(phy, phy_mat),
                                                     if(aces) list(aces_A_mat) else NULL),
             data = k,
             phy_match = phy_match,
             family = family,
             rate_model = rate_model,
             fit = fit, aces = aces,
             hyper = hyper,
             obs_error = obs_error,
             ...))
    })
    return(fits)
  } else {
    dat <- as.data.frame(dat)
  }

  namey <- data$node_name

  dat <- dat %>%
    dplyr::mutate(`root` = 1)


  tip_indexes <- 1:nrow(dat)
  node_indexes <- 1:ncol(phy_mat)

  A_mat <- phy_mat[namey, ]

  resp <- all.vars(formula[[2]])

  if(length(all.vars(formula)) > 1) {

    other_vars <- setdiff(all.vars(formula), resp)

    phy_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(dat))),
                                A = list(A_mat, 1),
                                effects = list(node_id = node_indexes,
                                               root = dat$root),
                                tag = "tces")

    other_stack <- INLA::inla.stack(data = list(y = dat[ , resp]),
                                A = list(A_mat, 1, 1),
                                effects = list(node_id = node_indexes,
                                               root = dat$root,
                                               covars = dat[ , other_vars, drop = FALSE]),
                                tag = "rates")

    phy_stack <- INLA::inla.stack(other_stack, phy_stack)
  } else {

    other_vars <- ""

    phy_stack <- INLA::inla.stack(data = list(y = dat[ , resp]),
                                A = list(A_mat, 1),
                                effects = list(node_id = node_indexes,
                                               root = dat$root),
                                tag = "rates")
  }

  if(aces) {

    aces_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(aces_A_mat))),
                                   A = list(aces_A_mat, 1),
                                   effects = list(node_id = 1:ncol(aces_A_mat),
                                                  root = rep(1, nrow(aces_A_mat))),
                                   tag = "aces")

    full_stack <- INLA::inla.stack(phy_stack, aces_stack)
  } else {
    full_stack <- phy_stack
    rm(phy_stack)
  }

  if(is.null(hyper)) {
    hyper <- "pc"
  }


  if(is.numeric(hyper)) {
    prior <- list(prec = list(initial = hyper, fixed = TRUE))
  } else {
    if(is.character(hyper)) {
      if(hyper == "pc") {
        dat_sd <- sd(dat[ , resp])
        l <- A_mat > 0
        e_var <- (sqrt(dat_sd / (sum(A_mat[l]) / length(phy$tip.label)) /
                         (sum(l) / length(phy$tip.label)))) * 3
        message("Automatically choosing prior for rate variance: exponential with 1% of probability density above ", e_var)
        prior <- list(prec = list(prior = "pc.prec", param = c(e_var, 0.01)))
      }
    }
    if(is.list(hyper)) {
      prior <- hyper
    }
  }

  obs_prior <- prior
  if(hyper == "pc") {
    obs_prior <- list(prec = list(prior = "pc.prec", param = c(e_var, 0.01)))
  }
  fam_cont <- switch(obs_error,
                     est = list(hyper = obs_prior),
                     one = list(hyper = list(prec = list(prior = "gaussian",
                                                         initial = 1,
                                                         fixed = TRUE))),
                     zero = list(hyper = list(prec = list(prior = "gaussian",
                                                          initial = 10,
                                                          fixed = TRUE))))


  if(rate_model == "ridge") {
    inla_form <- y ~ 0 + root + f(node_id, model = "iid",
                                  constr = FALSE,
                                  hyper = prior)


  }

  if(other_vars != "") {
    inla_form <- reformulate(c(deparse(formula[[3]]),
                               deparse(inla_form[[3]])),
                             "y")
  }

  message("Fitting model...")

  if(aces) {
    fit_modes <- INLA::inla(inla_form,
                       data = INLA::inla.stack.data(phy_stack),
                       family = family,
                       control.family = fam_cont,
                       control.predictor = list(A = INLA::inla.stack.A(phy_stack),
                                                compute = FALSE),
                       ...)

    fit <- INLA::inla(inla_form,
                      data = INLA::inla.stack.data(full_stack),
                      family = family,
                      control.family = fam_cont,
                      control.predictor = list(A = INLA::inla.stack.A(full_stack),
                                               compute = TRUE),
                      control.mode = list(theta = fit_modes$mode$theta, restart = FALSE),
                      ...)
  } else {
    fit <- INLA::inla(inla_form,
                      data = INLA::inla.stack.data(full_stack),
                      family = family,
                      control.predictor = list(A = INLA::inla.stack.A(full_stack),
                                               compute = TRUE))
  }

  nam <- rownames(fit$summary.fitted.values)
  rate_inds <- grep(".Predictor.", nam, fixed = TRUE)
  root_ind <- which(!is.na(full_stack$effects$data[ , "root"]))

  rate_index <- rate_inds[-root_ind]

  node_pred_index <- grep(".APredictor.", nam, fixed = TRUE)
  ace_ind <- INLA::inla.stack.index(full_stack, "aces")$data
  tip_ind <- setdiff(node_pred_index, ace_ind)
  tce_ind <- INLA::inla.stack.index(full_stack, "tces")$data

  attr(fit, "stack") <- full_stack
  attr(fit, "indexes") <- list(rates = rate_index,
                               node_predictions = node_pred_index,
                               aces = ace_ind,
                               tips = tip_ind,
                               tces = tce_ind)
  attr(fit, "node_names") <- list(tips = namey,
                                  nodes = if(aces) rownames(aces_A_mat) else NULL)

  fit

}
