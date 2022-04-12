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
#' @param family Family for the distribution of tip errors. Default is "gaussian", but can be any
#' supported by \code{inla}. See \code{names(inla.models()$likelihood)} for a list of possibilities
#' and use \code{\link[INLA]{inla.doc}} for details of all available families. Note that for
#' most families besides `"gaussian"`, the 'trait' who's evolution is being modeled will
#' be a 'latent' trait representing the 'mean' parameter of the error distribution and it will be
#' linked to the data through a standard transformation. As an example, a `"binomial"` model
#' with binary response will have a latent gaussian trait that represents the probability of the tip
#' traits being a one as opposed to a zero and will be transformed using a logit function (by default)
#' to keep the probability between 0 and 1.
#' @param rate_model The model to use for the evolutionary rate variation along the phylogeny
#' (in other words, the choice of prior for rates). Current choices are:
#' \itemize{
#' \item{`"iid"`}{Rates are completely independent along all branches but are shrunk
#' towards zero by an independent Normal prior. This is classic Bayesian Ridge Regression and has
#' a single hyperparameter whose prior determines the degree of shrinkage.}
#' \item{`"phylogenetic"`}{Rates are modeled using a covariance structure derived from the phylogeny.
#' Therefore, rates can be thought of as evolving on the phylogeny according to a Brownian motion
#' model. This limits the degree to which rates can differ from the rates on their parent branch.
#' This model can be helpful when trying to predict missing values at tree tips because it provides
#' more information on how to draw rate values for terminal branches. However, it may limit
#' the ability to detect rare rate shifts. In this case, one could consider a combination
#' of both `"phylogenetic"` and `"iid"` models.}
#' }
#' More rate models will be implemented in the future. Note that rate models can be combined
#' by passing a character vector of length greater than one to this argument (e.g.
#' c("iid", "phylogenetic")) will include both rate models as additive factors in the model.
#' @param rate_order Can be either `"first_order"` or `"second_order"`. If `"first_order"` rates
#' are the signed rate of evolution along each branch. If `"second_order"` rates are actually
#' the rate of rate change, that is we model rates of evolution along a branch as the rate
#' of the parent branch plus some deviation from this rate: the parameters of the model are these
#' rate deviations instead of the rates themselves.
#' @param fit
#' @param aces
#' @param hyper
#' @param obs_error
#' @param root
#' @param verbose Whether to print information about the progress of `fiber`
#' @param inla_verbose Whether [INLA::inla()] should be run in verbose mode.
#' @param threads Number of threads to use to parallelize computation of the root to
#' matrix (default: 1). This can speed things up quite a bit for large trees, but
#' also tends to use a lot of memory, because the entire tree must be copied to
#' each thread.
#' @param inla_threads How many threads should the inla fitting process use?
#' Default is `INLA::inla.getOption("num.threads")`. See [INLA::inla()] for how to set this correctly.
#' @param ... Additional arguments to be passed to [INLA::inla()]
#'
#' @return
#' @export
#'
#' @examples
fibre <- function(formula = ~ 1, phy, data = NULL,
                  phy_match = "auto",
                  family = "gaussian",
                  rate_ident = TRUE,
                  rate_comb = NULL,
                  rate_model = "iid",
                  rate_order = "first",
                  weights = NULL,
                  rate_weights = NULL,
                  aces = TRUE,
                  constr = TRUE,
                  hyper = NULL,
                  fit = TRUE,
                  obs_error = "est",
                  root = "est",
                  fixed_anc = NULL,
                  verbose = TRUE,
                  inla_verbose = FALSE,
                  threads = 1,
                  inla_threads = NULL,
                  ...) {

  assert_inla()

  if(is.null(inla_threads)) {
    inla_threads <- INLA::inla.getOption("num.threads")
  }

  rate_model <- match.arg(rate_model,
                          c("iid", "phylogenetic", "fixed"),
                          several.ok = TRUE)

  rate_order <- match.arg(rate_order,
                          c("first", "second"),
                          several.ok = TRUE)

  if(is.numeric(obs_error)) {
    fam_cont <- list(hyper = list(prec = list(prior = "gaussian",
                                                         initial = obs_error,
                                                         fixed = TRUE)))
  } else {
    if(is.character(obs_error) && obs_error == "est") {
      fam_cont <- list()
    } else {
      stop('obs_error must be a numeric constant or "est"" (for estimate)')
    }
  }

  # if(is.numeric(root)) {
  #   fam_cont <- list()
  # } else {
  #   if(is.character(obs_error) && obs_error == "est") {
  #     fam_cont <- list(hyper = list(prec = list(prior = "gaussian",
  #                                                        initial = obs_error,
  #                                                        fixed = TRUE)))
  #   } else {
  #     stop('obs_error must be a numeric constant or "est"" (for estimate)')
  #   }
  # }

  if(verbose) {
    message("Assembling model data and structure...")
  }


  if(inherits(phy, "phylo")) {
    message("Generating root-to-tip matrix...")
    if(aces) {
      phy_mat <- make_root2tip(phy, return_nodes = "both",
                               return_type = "list",
                               sparse = TRUE,
                               order = rate_order,
                               return_ages = "temporal" %in% rate_model,
                               threads = threads)
      if("temporal" %in% rate_model) {
        ages <- attr(phy_mat, "ages")
      }
      aces_A_mat <- phy_mat[[2]]
      phy_mat <- phy_mat[[1]]
    } else {
      phy_mat <- make_root2tip(phy, return_nodes = "tips",
                               sparse = TRUE,
                               order = rate_order,
                               return_ages = "temporal" %in% rate_model,
                               threads = threads)
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
      suppressMessages(fibre(formula = ~ 1, phy = c(list(phy, phy_mat),
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
    dplyr::mutate(`(Intercept/root)` = 1)

  tip_indexes <- 1:nrow(dat)
  node_indexes <- 1:ncol(phy_mat)

  A_mat <- phy_mat[namey, ]

  resp <- all.vars(formula[[2]])

  if(length(all.vars(formula)) > 1) {

    other_vars <- setdiff(all.vars(formula), resp)

    phy_stack <- INLA::inla.stack(data = list(y = dat[ , resp]),
                                A = list(A_mat, 1, 1),
                                effects = list(node_id = node_indexes,
                                               root = dat$`(Intercept/root)`,
                                               covars = dat[ , other_vars, drop = FALSE]),
                                tag = "rates")

  } else {

    other_vars <- ""

    phy_stack <- INLA::inla.stack(data = list(y = dat[ , resp]),
                                A = list(A_mat, 1),
                                effects = list(node_id = node_indexes,
                                               root = dat$`(Intercept/root)`),
                                tag = "rates")
  }

  if(aces) {

    tces_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(dat))),
                                A = list(A_mat, 1),
                                effects = list(node_id = node_indexes,
                                               root = dat$`(Intercept/root)`),
                                tag = "tces")

    aces_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(aces_A_mat))),
                                   A = list(aces_A_mat, 1),
                                   effects = list(node_id = 1:ncol(aces_A_mat),
                                                  root = rep(1, nrow(aces_A_mat))),
                                   tag = "aces")

    full_stack <- INLA::inla.stack(phy_stack, aces_stack, tces_stack)
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


  if(rate_model == "iid") {
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
                       verbose = inla_verbose,
                       ...)

    fit <- INLA::inla(inla_form,
                      data = INLA::inla.stack.data(full_stack),
                      family = family,
                      control.family = fam_cont,
                      control.predictor = list(A = INLA::inla.stack.A(full_stack),
                                               compute = TRUE),
                      control.mode = list(theta = fit_modes$mode$theta, restart = FALSE),
                      verbose = inla_verbose,
                      ...)
  } else {
    fit <- INLA::inla(inla_form,
                      data = INLA::inla.stack.data(full_stack),
                      family = family,
                      control.predictor = list(A = INLA::inla.stack.A(full_stack),
                                               compute = TRUE),
                      verbose = inla_verbose)
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

  class(fit) <- c("fibre_results", class(fit))

  fit

}


#' Phylogenetic Branch Random Effect. To be used in the `formula` argument to `fibre`.
#'
#' @param id Character or integer vector with ids to match data to phylogeny. If character,
#' elements are matched to tip labels in `phy`. If integer, they will treated as indexes into
#' `phy$ip.label`. Can also be a [ape::phylo] object, in which case the `phy` argument is
#' ignored, and the `id` will be an integer indexing the tip labels (in other words, this
#' works if you know your data is already in the same order as your phylogeny's tip labels).
#' @param phy A [ape::phylo] object with the phylogeny to be used in the `phybre` model.
#' @param rate_model
#' @param rate_order
#' @param weights
#' @param rate_weights
#' @param predict_effects
#' @param constr
#' @param hyper
#' @param multiple_order
#' @param fixed_anc
#'
#' @return
#' @export
#'
#' @examples
p <- function(id,
              phy = NULL,
              rate_model = c("iid", "phylogenetic", "fixed"),
              rate_order = c("first", "second"),
              weights = NULL,
              rate_weights = NULL,
              aces = TRUE,
              constr = TRUE,
              hyper = NULL,
              fixed_anc = NULL,
              rate_group = 1L,
              multivariate = 0,
              name = NULL) {


  rate_model <- match.arg(rate_model)
  rate_order <- match.arg(rate_order)

  if(inherits(id, "phylo")) {
    phy <- id
    id <- seq_along(id$tip.label)
  }

  if(rate_model == "fixed" && is.null(rate_weights)) {
    stop("rate_model of 'fixed' requires rate_weights to be specified")
  }

  rate_group <- as.integer(rate_group)

  if(is.null(name)) {
    name <- make_name(rate_model, rate_order, weights, rate_weights)
  }

  multiple_order <- .fibre_env$multiple_order

  if(!is.null(rate_weights) && rate_weights == "age") {
    age <- TRUE
  } else {
    age <- FALSE
  }

  if(multiple_order) {
    rtp_mat <- make_root2tip(phy,
                             return_nodes = "both",
                             return_type = "list",
                             order = "both",
                             return_ages = age
                             )
    parents <- attr(rtp_mat, "parents")
    if(age) {
      ages <- attr(rtp_mat, "ages")
    }
    rtp_mat <- rtp_mat[[rate_order]]
  } else {
    rtp_mat <- make_root2tip(phy,
                             return_nodes = "both",
                             return_type = "list",
                             order = rate_order,
                             return_ages = age
                             )
    parents <- attr(rtp_mat, "parents")
    if(age) {
      ages <- attr(rtp_mat, "ages")
    }
  }

  if(age) {
    rate_weights <- ages
  }

  if(!is.null(weights)) {
    rtp_mat <- lapply(rtp_mat, function(x) x * weights)
  }

  if(rate_model != "fixed" && !is.null(rate_weights)) {
    rtp_mat <- lapply(rtp_mat, function(x) t(t(x) * rate_weights))
  }

  if(aces) {
    rtp_mat_tces <- rtp_mat$tips
    rtp_mat_aces <- rtp_mat$internal

    tces_effect <- list(node_id = 1:ncol(rtp_mat_tces))
    aces_effect <- list(node_id = 1:ncol(rtp_mat_aces))

    names(tces_effect) <- name
    names(aces_effect) <- name

    tces_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(rtp_mat_tces))),
                                   A = list(rtp_mat_tces),
                                   effects = tces_effect,
                                   tag = "tces")

    aces_stack <- INLA::inla.stack(data = list(y = rep(NA, nrow(rtp_mat_aces))),
                                   A = list(rtp_mat_aces),
                                   effects = aces_effect,
                                   tag = "aces")

  } else {
    tces_stack <- NULL
    aces_stack <- NULL
  }

  rtp_mat <- rtp_mat$tips
  rtp_mat <- rtp_mat[id, ]

  if(rate_model == "phylogenetic") {
    Cmat <- make_Cmatrix(phy)
  } else {
    Cmat <- NULL
  }

  if(constr) {

    extra_constr <- make_extraconstr(phy, parents, rate_order)

  } else {
    extra_constr <- list(A = NULL, e = NULL)
  }

  if(rate_model == "phylogenetic") {
    rate_model <- "generic0"
  }

  if(rate_model == "fixed") {
    data <- list(fixed_effect = rate_weights)
  } else {
    data <- list(node_id = 1:ncol(rtp_mat))
  }

  names(data) <- name

  list(rtp_mat = rtp_mat,
       tces_stack = tces_stack,
       aces_stack = aces_stack,
       Cmat = Cmat,
       extra_constr = extra_constr,
       data = data,
       hyper = hyper,
       rate_model = rate_model,
       name = gsub("phy_", "", name))

}
