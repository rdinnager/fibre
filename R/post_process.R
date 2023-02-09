fibre_rates <- function(fibre_mod, return_sample = 0, return_summary = return_sample == 0,
                        return_marginal = return_sample == 0, 
                        return_type = c("pf", "data.frame")) {
  rates <- fibre_mod$random
  
  if(return_sample > 0) {
    if(fibre_mod$engine != "inla") {
      rlang::abort("Returning a sample from the posterior is only supported when engine = 'inla'")
    }
    samps <- lapply(rates, function(x) purrr::map(x$marginal,
                                                  ~ INLA::inla.rmarginal(return_sample, .x)) %>%
                      do.call(rbind, .))
    
    samps <- dplyr::as_tibble(samps)
    colnames(samps) <- paste0("rate_sample_", colnames(samps))
    rates <- rates %>%
      dplyr::bind_cols(sample = samps)
  }
  if(!return_marginal) {
    if(fibre_mod$engine == "inla") {
      rates$marginal <- NULL
    }
  }
}