#' @export
fibre_formula_blueprint <- function(intercept = FALSE, allow_novel_levels = FALSE, 
                                    indicators = "traditional", 
                                    composition = "tibble",
                                    fixed_blueprint = NULL) {
  
  new_fibre_formula_blueprint(intercept = intercept, 
                              allow_novel_levels = allow_novel_levels, 
                              indicators = indicators, 
                              composition = composition,
                              fixed_blueprint = fixed_blueprint)  
  
}

new_fibre_formula_blueprint <- function(intercept = FALSE, allow_novel_levels = FALSE, 
                                        indicators = "traditional",
                                        composition = "tibble",
                                        fixed_blueprint = NULL,
                                        ...) {
  
  hardhat::new_default_formula_blueprint(intercept = intercept, 
                                         allow_novel_levels = allow_novel_levels, 
                                         indicators = indicators, 
                                         composition = composition,
                                         fixed_blueprint = fixed_blueprint,
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
  
  blueprint <- hardhat::update_blueprint(blueprint, fixed_blueprint = dat$blueprint)
  
  info$new_form <- NULL
  dat$extras$model_info <- info$bre_info
  list(predictors = dat$predictors,
       outcomes = dat$outcomes,
       blueprint = blueprint,
       extras = dat$extras)
  
}

#' @importFrom hardhat run_forge
#' @export
run_forge.fibre_formula_blueprint <- function(blueprint, new_data, ..., outcomes = FALSE) {
  
  form <- blueprint$formula  
  info <- parse_formula_for_mold(form, data = new_data)
  #blueprint_default <- do.call(hardhat::new_default_formula_blueprint, as.list(blueprint$fixed_blueprint))
  
  new_dat <- hardhat::forge(new_data = new_data,
                            blueprint = blueprint$fixed_blueprint,
                            outcomes = outcomes)
  info$new_form <- NULL
  new_dat$extras$model_info <- info$bre_info
  list(predictors = new_dat$predictors,
       outcomes = new_dat$outcomes,
       extras = new_dat$extras)
  
}

#' @importFrom hardhat refresh_blueprint
#' @export
refresh_blueprint.fibre_formula_blueprint <- function(blueprint) {
  do.call(new_fibre_formula_blueprint, as.list(blueprint))
}


