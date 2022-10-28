#' Create a Evolutionary Autodecoder Model
#'
#' @param latent_dim 
#' @param decoder A `torch::nn_model()` specifying a 'decoder' network architecture. 
#' The decoder network should accept a 2 dimensional `torch::torch_tensor()` with
#' first d
#'
#' @return
#' @export
#'
#' @examples
evo_autodecoder <- function(latent_dim, n_edges, decoder, reconstruction_loss, device, decoder_args = list(), loss_args = list()) {
  torch::nn_module(
    "EAD",

    initialize = function(latent_dim, n_edges) {
        self$latent_dim <- latent_dim
        self$n_edges <- n_edges
        
        self$latent_rate_means <- torch::nn_parameter(torch::nn_init_normal_(torch::torch_empty(n_edges,
                                                                                                latent_dim)))
        
        self$latent_rate_log_vars <- torch::nn_parameter(torch::nn_init_normal_(torch::torch_empty(latent_dim, 
                                                                                                n_edges)))
        self$encoder <- function(x, rates) torch::nnf_linear(x, rates)

        self$decoder <- decoder
    },

    decode = function(z, ...) {
        self$decoder(z, ...)
    },
    
    encode = function(x, rates) {
      self$encode(x, rates)
    },

    reparameterize = function(mean, log_var) {
        std <- torch_tensor(0.5, device = device) * log_var
        eps <- torch_randn_like(std)
        eps * std + mean
    },

    loss_function = function(reconstruction, observed) {
        reconstruction_loss <- rlang::exec(self$reconstruction_loss, reconstruction = reconstruction, observed = observed, !!!loss_args)
        kl_loss <- torch_tensor(-0.5, device = device) * torch_sum(torch_tensor(1, device = "cuda") + self$latent_rates_log_vars - self$latent_rates_means^2 - self$latent_rates_log_vars$exp())
        loss <- reconstruction_loss + kl_loss
        list(loss, reconstruction_loss, kl_loss)
    },

    forward = function(x) {
        rates <- self$reparameterize(self$latent_rate_means, 
                                     self$latent_rates_log_vars)
        z <- self$encode(x, rates)
        recon <- rlang::exec(self$decode, z = z, !!!decoder_args)
        list(recon, x)
    },

    sample = function(num_samples, current_device) {
        z <- torch_randn(num_samples, self$latent_dim)
        z <- z$to(device = current_device)
        samples <- rlang::exec(self$decode, z = z, !!!decoder_args)
        samples
    }
  )  
  
}

