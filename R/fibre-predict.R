#' Predict from a `fibre`
#'
#' @param object A `fibre` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"numeric"` for numeric predictions.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @examples
#' train <- mtcars[1:20,]
#' test <- mtcars[21:32, -1]
#'
#' # Fit
#' mod <- fibre(mpg ~ cyl + log(drat), train)
#'
#' # Predict, with preprocessing
#' predict(mod, test)
#'
#' @export
predict.fibre <- function(object, new_data = NULL, type = "numeric", ...) {
  if(is.null(new_data)) {
    return(object$saved_predictions %>% dplyr::select(-dplyr::all_of("label")))
  }
  
  forged <- hardhat::forge(new_data, object$blueprint)
  pfcs <- purrr::map(forged$extras$model_info,
                     "phyf")
  to_pred <- purrr::map(pfcs,
                        phyf::pf_labels) %>%
    dplyr::bind_cols(.name_repair = ~ vctrs::vec_as_names(..., 
                                                          repair = "unique", 
                                                          quiet = TRUE)) %>%
    dplyr::bind_cols(forged$predictors, 
                     .name_repair = ~ vctrs::vec_as_names(..., 
                                                          repair = "unique", 
                                                          quiet = TRUE)) %>%
    tidyr::unite(label, dplyr::everything())
  
  preds <- to_pred %>%
    dplyr::left_join(object$saved_predictions, by = "label") %>%
    dplyr::select(-dplyr::all_of("label"))
  
  to_pred <- to_pred$label[!to_pred$label %in% object$saved_predictions$label]
  
  if(length(to_pred) == 0) {
    return(preds)
  }
  rlang::arg_match(type, valid_fibre_predict_types())
  predict_fibre_bridge(type, object, forged$predictors)
}

valid_fibre_predict_types <- function() {
  c("numeric")
}

# ------------------------------------------------------------------------------
# Bridge

predict_fibre_bridge <- function(type, model, predictors) {
  predictors <- as.matrix(predictors)

  predict_function <- get_fibre_predict_function(type)
  predictions <- predict_function(model, predictors)

  hardhat::validate_prediction_size(predictions, predictors)

  predictions
}

get_fibre_predict_function <- function(type) {
  switch(
    type,
    numeric = predict_fibre_numeric
  )
}

# ------------------------------------------------------------------------------
# Implementation

predict_fibre_numeric <- function(model, predictors) {
  predictions <- rep(1L, times = nrow(predictors))
  hardhat::spruce_numeric(predictions)
}
