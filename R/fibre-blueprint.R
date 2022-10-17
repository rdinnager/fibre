#' @export
fibre_formula_blueprint <- function(intercept = FALSE, allow_novel_levels = FALSE, 
                                    indicators = "traditional", 
                                    composition = "tibble") {
  
  new_fibre_formula_blueprint(intercept = intercept, 
                              allow_novel_levels = allow_novel_levels, 
                              indicators = indicators, 
                              composition = composition)  
  
}

new_fibre_formula_blueprint <- function(intercept = FALSE, allow_novel_levels = FALSE, 
                                        indicators = "traditional",
                                        composition = "tibble",
                                        ...) {
  
  hardhat::new_default_formula_blueprint(intercept = intercept, 
                                         allow_novel_levels = allow_novel_levels, 
                                         indicators = indicators, 
                                         composition = composition,
                                         subclass = "fibre_formula_blueprint",
                                         ...)
  
}

#' @importFrom hardhat run_mold
#' @export
run_mold.fibre_formula_blueprint <- function(blueprint, data, ...) {
  
  form <- blueprint$formula  
  info <- parse_formula_for_mold(form, data = data)
  blueprint_default <- do.call(hardhat::new_default_formula_blueprint, as.list(blueprint))
  dat <- hardhat::mold(info$new_form, data = data, 
                       blueprint = blueprint_default)
  info$new_form <- NULL
  dat$extras$model_info <- info$ps_info
  list(predictors = dat$predictors,
       outcomes = dat$outcomes,
       blueprint = dat$blueprint,
       extras = dat$extras)
  
}

#' @importFrom hardhat refresh_blueprint
#' @export
refresh_blueprint.fibre_formula_blueprint <- function(blueprint) {
  do.call(new_fibre_formula_blueprint, as.list(blueprint))
}


