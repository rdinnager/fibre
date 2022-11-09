parse_formula_for_mold <- function(form, data = NULL, debug = FALSE) {

  ts <- terms.formula(form, specials = c("bre",
                                         "bre_brownian",
                                         "bre_second_order",
                                         "re"),
                      keep.order = TRUE)
  
  bres <- rownames(attr(ts, "factors"))[unlist(attr(ts, "specials"))]
  bre_ind <- which(attr(ts, "term.labels") %in% bres)
  
  ts <- ts[-bre_ind]
  
  new_form <- formula(ts)
  
  if(as.character(new_form)[3] == "1") {
    new_form <- as.formula(rlang::expr(!!rlang::f_lhs(new_form) ~ NULL))
  }
  
  bre_info <- lapply(bres, function(x) rlang::eval_tidy(rlang::parse_expr(x),
                                                       data = data))

  list(new_form = new_form, bre_info = bre_info)

}

#' Specify a branch length (random) effect
#' 
#' This function is meant to be called only in the `formula` argument of `fibre()`.
#'
#' @param phyf A `pfc` column containing the phylogenetic structure
#' @param rate_distribution What distribution to use to model rates of evolution?
#' @param hyper Hyper parameters as a list. Specify the prior distribution for 
#' `engine = INLA` models here. Default is a penalised complexity prior with 10%
#' prior probability density greater than 1, which tend to work well for 
#' standardised Gaussian responses and Binomial responses.
#' @param latent How many latent variables to generate in `engine = INLA` models.
#' Default is none.
#' @param label An optional label used to identify the random effect later 
#' The default is a label generated from the expression in `phyf`
#' @param standardise Should the `pfc` object be standardised based on
#' it's implied typical variance for terminal nodes? Default: `TRUE`.
#' This helps random effects to be comparable to each other.
#'
#' @return A list of data to be used by the model.
#' @export
bre <- function(phyf, 
                rate_distribution = c("iid", "laplacian", "student-t", "horseshoe", "Brownian", "re"),
                hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))),
                latent = 0,
                label = NULL,
                standardise = TRUE) {
  
  if(is.null(label)) label <- rlang::expr_label(substitute(phyf))
  
  if(standardise) {
    standard <- mean(phyf::pf_flow_sum(phyf))
    phyf <- phyf / standard
  }
  
  rate_dist <- match.arg(rate_distribution)
  list(phyf = phyf,
       rate_dist = rate_dist,
       hyper = hyper,
       latent = latent,
       label = label,
       standard = standard)
  
}

#' Specify a random effect
#' 
#' This function is meant to be called only in the `formula` argument of `fibre()`.
#'
#' @param groups A character or factor column containing the grouping variable 
#' for the random effect
#' @inheritParams bre
#'
#' @return A list of data to be used by the model.
#' @export
re <- function(groups, 
               hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))),
               label = NULL,
               standardise = TRUE) {
  
  if(is.null(label)) label <- rlang::expr_label(substitute(groups))
  phyf <- pf_as_pfc(groups)
  
  bre(phyf = phyf,
      rate_dist = "re",
      hyper = hyper,
      latent = 0,
      label = label,
      standardise = TRUE)
  
}

#' Specify a branch length (random) effect for a Brownian motion model
#' 
#' This function is meant to be called only in the `formula` argument of `fibre()`.
#'
#' @inheritParams bre
#'
#' @return A list of data to be used by the model.
#' @export
bre_brownian <- function(phyf, 
                         hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))),
                         latent = 0,
                         label = NULL,
                         standardise = TRUE) {
  
  if(is.null(label)) label <- rlang::expr_label(substitute(phyf))
  
  bre(phyf = sqrt(phyf),
      rate_dist = "Brownian",
      hyper = hyper,
      latent = latent,
      label = label,
      standardise = TRUE)
  
}

#' Specify a branch length (random) effect for a 'Second Order' Brownian 
#' motion model
#' 
#' This function is meant to be called only in the `formula` argument of `fibre()`.
#'
#' @inheritParams bre
#'
#' @return A list of data to be used by the model.
#' @export
bre_second_order <- function(phyf, 
                             hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))),
                             latent = 0,
                             label = NULL,
                             standardise = TRUE) {
  
  if(is.null(label)) label <- rlang::expr_label(substitute(phyf))
  
  bre(phyf = sqrt(phyf::pf_second_order(phyf)),
      rate_dist = "iid",
      hyper = hyper,
      latent = latent,
      label = label,
      standardise = TRUE)
  
}


make_inla_formula <- function(dat, y) {
  
  preds <- colnames(dat)
  re <- grep("_indexes", preds)
  latent <- grep("latent_", preds)
  copies <- grep("copy_latent_", preds)
  
  if(length(latent) > 0) {
    re <- setdiff(re, latent)
    latent = setdiff(latent, copies)
  }
  
  #preds <- purrr::map_chr(rlang::syms(preds), rlang::expr_label)
  
  ys <- colnames(y)
  #ys <- purrr::map_chr(rlang::syms(ys), rlang::expr_label)
  
  if(length(re) > 0) {
    hypers_re <- paste0("hyper_re_", seq_along(re))
    re_fs <- glue::glue("f({preds[re]}, model = 'iid', hyper = {hypers_re})")
  } else {
    hypers_re <- character(0)
  }
  
  if(length(latent) > 0) {
    hypers_latent <- paste0("hyper_latent_", seq_along(latent))
    latent_fs <- glue::glue("f({preds[latent]}, model = 'iid', hyper = {hypers_latent})")
  }
  
  if(length(latent) > 0) {
    hypers_copy <- paste0("hyper_copy_", seq_along(copies))
    copies_of <- sapply(strsplit(preds[copies], ".", fixed = TRUE), function(x) x[2])
    copies_of <- gsub("copy_", "", copies_of)
    copy_fs <- glue::glue('f({preds[copies]}, copy = "{copies_of}", hyper = {hypers_copy})')
  } else {
    hypers_latent <- character(0)
    hypers_copy <- character(0)
  }
  
  list(form = as.formula(paste(
    ifelse(length(ys) > 1, "y", paste(ys, collapse = " + ")),
    "~ -1 +",
    paste(preds[-c(re, latent, copies)], collapse = " + "),
    ifelse(length(re) > 0, paste0(" + ", paste(re_fs, collapse = " + ")), ""),
    ifelse(length(latent) > 0, paste0(c(" + "),
                                     paste(latent_fs, collapse = " + "),
                                     " + ",
                                     paste(copy_fs, collapse = " + ")), 
           "")
    )
  ),
  hypers = list(re = hypers_re, latent = hypers_latent, copy = hypers_copy))
}


