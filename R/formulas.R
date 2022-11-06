parse_formula_for_mold <- function(form, data = NULL, debug = FALSE) {

  ts <- terms.formula(form, specials = c("bre",
                                         "bre_brownian",
                                         "bre_second_order"),
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
#' `engine = INLA` models here.
#' @param latent How many latent variables to generate in `engine = INLA` models.
#' Default is none.
#'
#' @return A list of data to be used by the model.
#' @export
bre <- function(phyf, 
                rate_distribution = c("iid", "laplacian", "student-t", "horseshoe", "Brownian"),
                hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))),
                latent = 0) {
  
  rate_dist <- match.arg(rate_distribution)
  list(phyf = phyf,
       rate_dist = rate_dist,
       hyper = hyper,
       latent = latent)
  
}

#' Specify a branch length (random) effect for a Brownian motion model
#' 
#' This function is meant to be called only in the `formula` argument of `fibre()`.
#'
#' @param phyf A `pfc` column containing the phylogenetic structure
#' @param hyper Hyper parameters as a list. Specify the prior distribution for 
#' `engine = INLA` models here.
#' @param latent How many latent variables to generate in `engine = INLA` models.
#' Default is none.
#'
#' @return A list of data to be used by the model.
#' @export
bre_brownian <- function(phyf, 
                         hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))),
                         latent = 0) {
  
  bre(phyf = sqrt(phyf),
      rate_dist = "Brownian",
      hyper = hyper,
      latent = latent)
  
}

#' Specify a branch length (random) effect for a 'Second Order' Brownian 
#' motion model
#' 
#' This function is meant to be called only in the `formula` argument of `fibre()`.
#'
#' @param phyf A `pfc` column containing the phylogenetic structure
#' @param hyper Hyper parameters as a list. Specify the prior distribution for 
#' `engine = INLA` models here.
#' @param latent How many latent variables to generate in `engine = INLA` models.
#' Default is none.
#'
#' @return A list of data to be used by the model.
#' @export
bre_second_order <- function(phyf, 
                             hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))),
                             latent = 0) {
  
  bre(phyf = sqrt(phyf::pf_second_order(phyf)),
      rate_dist = "iid",
      hyper = hyper,
      latent = latent)
  
}


make_inla_formula <- function(dat, y) {
  
  preds <- colnames(dat)
  re <- grep("_indexes", preds)
  hypers <- paste0("hyper_", seq_along(re))
  #preds <- purrr::map_chr(rlang::syms(preds), rlang::expr_label)
  
  ys <- colnames(y)
  #ys <- purrr::map_chr(rlang::syms(ys), rlang::expr_label)
  
  fs <- glue::glue("f({preds[re]}, model = 'iid', hyper = {hypers})")
  
  as.formula(paste(
    paste(ys, collapse = " + "),
    "~ -1 +",
    paste(preds[-re], collapse = " + "),
    "+",
    paste(fs, collapse = " + ")
    )
  )
}


