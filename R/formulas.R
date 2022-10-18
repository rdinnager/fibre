parse_formula_for_mold <- function(form, data = NULL, debug = FALSE) {

  ts <- terms.formula(form, specials = "bre", keep.order = TRUE)
  
  bres <- rownames(attr(ts, "factors"))[attr(ts, "specials")$bre]
  bre_ind <- which(attr(ts, "term.labels") %in% bres)
  
  ts <- ts[-bre_ind]
  
  new_form <- formula(ts)
  
  bre_info <- lapply(bres, function(x) rlang::eval_tidy(rlang::parse_expr(x),
                                                       data = data))

  list(new_form = new_form, bre_info = bre_info)

}

bre <- function(phyf, 
                rate_distribution = c("iid", "gaussian", "ridge", "laplacian", "lasso", "double-exponential", "student-t", "horseshoe", "diverging", "fixed", "Brownian"),
                hyper = NULL,
                latent = 0) {
  
  rate_dist <- match.arg(rate_distribution)
  list(phyf = phyf,
       rate_dist = rate_dist,
       hyper = hyper,
       latent = latent)
  
}