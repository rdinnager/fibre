new_fibre <- function(fixed, random, hyper, model, saved_predictions, engine, blueprint) {
  ob <- hardhat::new_model(fixed = fixed,
                     random = random,
                     hyper = hyper,
                     model = model,
                     saved_predictions = saved_predictions,
                     engine = engine,
                     blueprint = blueprint, 
                     class = "fibre")
}

#' @export
print.fibre <- function(x, n = 10, ...) {
  
  cli::cli_text("A <{class(x)[1]}> model with {length(x$random)} phylogenetic random {cli::qty(length(x$random))} effect{?s}")
  
  cat("\n")
  
  cli::cli_rule(left = cli::style_bold("Fixed Effects"))
  
  cat("\n")
  
  if(x$engine == "inla") {
    print(x$fixed %>%
          dplyr::mutate(marginal = spark_hist_with_padding(marginal)))
  }
  if(x$engine == "glmnet") {
    print(x$fixed %>%
          dplyr::mutate(dotplot = spark_dotplot(coef)))
  }
  
  cat("\n")
  
  cli::cli_rule(left = cli::style_bold("Random Effect Hyper-Parameters"))
  
  cat("\n")
  
  if(x$engine == "inla") {
    print(x$hyper %>%
            dplyr::mutate(marginal = spark_hist_with_padding(marginal)))
  }
  if(x$engine == "glmnet") {
    print(x$hyper %>%
            dplyr::mutate(dotplot = spark_dotplot(value)))
  }
  
  cat("\n")
  
  for(i in seq_along(x$random)) {
    
    cli::cli_rule(left = cli::style_bold("Random Effects ('Rates') for {names(x$random)[i]}"))
    cli::cli_text("A Random Effect with {nrow(x$random[[i]])} estimated rates.")
    
    cat("\n")
    
    cli::cli_text("First {n} highest rates: ")
    
    if(x$engine == "inla") {
      print(x$random[[i]] %>%
              dplyr::slice_max(abs(mean), n = n) %>%
              dplyr::mutate(marginal = spark_hist_with_padding(marginal)))
    }
    if(x$engine == "glmnet") {
      print(x$random[[i]] %>%
              dplyr::slice_max(abs(coef), n = n) %>%
              dplyr::mutate(dotplot = spark_dotplot(coef)))
    }
    cli::cli_text("... with {nrow(x$random[[i]]) - n} more rate{?s}.")
    cli::cli_alert_info("Use `print(n = ...)` to see more rates.")
    cat("\n")
  }
}