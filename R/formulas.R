parse_formula_for_mold <- function(form, data = NULL, debug = FALSE) {

  ts <- terms.formula(form, specials = "p", keep.order = TRUE)
  
  ps <- rownames(attr(ts, "factors"))[attr(ts, "specials")$p]
  ps_ind <- which(attr(ts, "term.labels") %in% ps)
  
  ts <- ts[-ps_ind]
  
  new_form <- formula(ts)
  
  ps_info <- lapply(ps, function(x) rlang::eval_tidy(rlang::parse_expr(x),
                                                     data = data))

  list(new_form = new_form, ps_info = ps_info)

}

p <- function(phyf, 
              rate_distribution = c("iid", "gaussian", "ridge", "laplacian", "lasso", "double-exponential", "student-t", "horseshoe", "diverging", "fixed", "Brownian"),
              hyper = NULL,
              latent = 0,
              mixture_of = NULL) {
  
  rate_dist <- match.arg(rate_distribution)
  list(phyf = phyf,
       rate_dist = rate_dist,
       hyper = hyper,
       latent = latent,
       mixture_of = mixture_of)
  
}