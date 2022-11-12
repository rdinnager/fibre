#' Load a model
#'
#' @param name Name of the model. Currently only `"bird_beak"`.
#'
#' @return A `torch::nn_module()` with pre-trained weights
#' @export
#'
#' @examples
#' if(torch::torch_is_installed) {
#' model <- load_model("bird_beaks")
#' }
load_model <- function(name) {
  
  switch(name,
         bird_beaks = load_bird_beak_model())
  
  
}

load_bird_beak_model <- function() {
  sd1 <- torch::load_state_dict(system.file("models/sdf_net.pt", package = "fibre"))
  sdfnet <- sdf_net()
  sdfnet$load_state_dict(sd1)
  sdfnet$eval()
  sdfnet
}