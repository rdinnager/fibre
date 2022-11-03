new_fibre <- function(fixed, random, hyper, model, saved_predictions, blueprint) {
  hardhat::new_model(fixed = fixed,
                     random = random,
                     hyper = hyper,
                     model = model,
                     saved_predictions = saved_predictions,
                     blueprint = blueprint, 
                     class = "fibre")
}
