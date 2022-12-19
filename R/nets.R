#' A signed distance field based neural network model for generating 3d shapes
#'
#' @param n_latent Number of dimensions for the latent space
#' @param breadth Breadth of the multilayer perceptron networks
#'
#' @return A `torch::nn_module()`
#' @export
#'
#' @examples
#' if(torch::torch_is_installed()) {
#' sdf_net()
#' }
sdf_net <- torch::nn_module("sdf_net",

                            initialize = function(n_latent = 64, breadth = 512) {

                              self$layers1 <- torch::nn_sequential(
                                torch::nn_linear(in_features = 3L + n_latent,
                                                 out_features = breadth),
                                torch::nn_relu(inplace = TRUE),
                                torch::nn_linear(in_features = breadth,
                                                 out_features = breadth),
                                torch::nn_relu(inplace = TRUE),
                                torch::nn_linear(in_features = breadth,
                                                 out_features = breadth),
                                torch::nn_relu(inplace = TRUE)
                              )

                              self$layers2 <- torch::nn_sequential(
                                torch::nn_linear(in_features = breadth + n_latent + 3,
                                                 out_features = breadth),
                                torch::nn_relu(inplace = TRUE),
                                torch::nn_linear(in_features = breadth,
                                                 out_features = breadth),
                                torch::nn_relu(inplace = TRUE),
                                torch::nn_linear(in_features = breadth,
                                                 out_features = 1L),
                                torch::nn_tanh()
                              )

                            },

                            forward = function(points, latent_codes) {
                              if(latent_codes$shape[1] == 1) {
                                latent_codes <- latent_codes$`repeat`(c(points$shape[1], 1L))
                              }
                              input <- torch::torch_cat(list(points, latent_codes), dim = 2L)
                              x <- self$layers1(input)
                              x <- torch::torch_cat(list(x, input), dim = 2L)
                              x <- self$layers2(x)

                              return(x$squeeze(2L))
                            },

                            get_activations = function(points, latent_codes) {
                              if(latent_codes$shape[1] == 1) {
                                latent_codes <- latent_codes$`repeat`(c(points$shape[1], 1L))
                              }
                              input <- torch::torch_cat(list(points, latent_codes), dim = 2L)
                              x <- self$layers1(input)
                              x <- torch::torch_cat(list(x, input), dim = 2L)
                              x <- self$layers2(x)

                              return(x$squeeze(2L))
                            },

                            get_normals = function(points, latent_code = NULL) {
                              # if(latent_code$requires_grad | points$requires_grad) {
                              #   stop('get_normals may only be called with tensors that don\'t require grad.')
                              # }
                              points$requires_grad <- TRUE
                              if(is.null(latent_code)) {
                                sdf <- self$forward(points)
                              } else {
                                if(latent_code$shape[1] == 1) {
                                  latent_code <- latent_code$`repeat`(c(points$shape[1], 1L))
                                }
                                #latent_code$requires_grad <- TRUE
                                sdf <- self$forward(points, latent_code)
                              }

                              ## needed to get torch to release GPU memory
                              ## this code hangs after a few iterations otherwise
                              ## but adds a lot of overhead
                              gc()

                              normals <- torch::autograd_grad(sdf, points,
                                                              torch::torch_ones_like(sdf))
                              normals <- normals[[1]] / torch::torch_norm(normals[[1]], dim = 2L, keepdim = TRUE)
                              return(normals)

                            },

                            get_normals_in_batches = function(points, latent_code = NULL,
                                                              batch_size = 10000,
                                                              return_cpu = TRUE,
                                                              cuda = FALSE) {
                              if(is.null(latent_code)) {
                                if(!cuda) {
                                  res <- map_cat(torch::torch_split(points, batch_size),
                                                   ~ self$get_normals(.x))
                                } else {
                                  if(return_cpu) {
                                    res <- map_cat(torch::torch_split(points, batch_size),
                                                   ~ self$get_normals(.x$cuda())$cpu())
                                  } else {
                                    res <- map_cat(torch::torch_split(points, batch_size),
                                                   ~ self$get_normals(.x$cuda()))
                                  }
                                }
                              } else {
                                if(!cuda) {
                                  res <- map_cat(torch::torch_split(points, batch_size),
                                                   ~ self$get_normals(.x, latent_code))
                                } else {
                                  if(return_cpu) {
                                    res <- map_cat(torch::torch_split(points, batch_size),
                                                   ~ self$get_normals(.x$cuda(), latent_code$cuda())$cpu())
                                  } else {
                                    res <- map_cat(torch::torch_split(points, batch_size),
                                                   ~ self$get_normals(.x$cuda(), latent_code$cuda()))
                                  }
                                }
                              }

                              res
                            },

                            evaluate_in_batches = function(points, latent_code = NULL,
                                                           batch_size = 10000, return_cpu = TRUE,
                                                           cuda = FALSE) {

                              if(is.null(latent_code)) {
                                if(!cuda) {
                                  torch::with_no_grad({
                                      res <- map_cat(torch::torch_split(points, batch_size),
                                                      ~ self$forward(.x))
                                    })
                                } else {
                                  if(return_cpu) {
                                    torch::with_no_grad({
                                         res <- map_cat(torch::torch_split(points, batch_size),
                                                        ~ self$forward(.x$cuda())$cpu())
                                    })
                                  } else {
                                    torch::with_no_grad({
                                      res <- map_cat(torch::torch_split(points, batch_size),
                                                     ~ self$forward(.x$cuda()))
                                    })
                                  }
                                }
                              } else {
                                if(!cuda) {
                                  torch::with_no_grad({
                                      res <- map_cat(torch::torch_split(points, batch_size),
                                                     ~ self$forward(.x, latent_code))
                                    })
                                } else {
                                  if(return_cpu) {
                                    torch::with_no_grad({
                                      res <- map_cat(torch::torch_split(points, batch_size),
                                                     ~ self$forward(.x$cuda(), latent_code$cuda())$cpu())
                                    })
                                  } else {
                                    torch::with_no_grad({
                                      res <- map_cat(torch::torch_split(points, batch_size),
                                                     ~ self$forward(.x$cuda(), latent_code))
                                    })
                                  }
                                }
                              }

                              res
                            },

                            render_image = function(latent_code = NULL,
                                                    resolution = 800,
                                                    camera_position = default_camera(),
                                                    light_position = default_light(),
                                                    threshold = 0.0005,
                                                    sdf_offset = 0,
                                                    iterations = 1000,
                                                    ssaa = 2,
                                                    radius = 1.0,
                                                    crop = FALSE,
                                                    color = c(R = 255 / 255, G = 237 / 255, B = 95 / 255),
                                                    vertical_cutoff = NULL,
                                                    max_ray_move = 0.05,
                                                    plot = TRUE,
                                                    cuda = FALSE,
                                                    batch_size = 50000,
                                                    verbose = FALSE) {

                                                            render_image(self, latent_code = latent_code,
                                                                         resolution = resolution,
                                                                         camera_position = camera_position,
                                                                         light_position = light_position,
                                                                         threshold = threshold,
                                                                         sdf_offset = sdf_offset,
                                                                         iterations = iterations,
                                                                         ssaa = ssaa,
                                                                         radius = radius,
                                                                         crop = crop,
                                                                         color = color,
                                                                         vertical_cutoff = vertical_cutoff,
                                                                         max_ray_move = max_ray_move,
                                                                         plot = plot,
                                                                         cuda = cuda,
                                                                         batch_size = batch_size,
                                                                         verbose = verbose)

                            },

                            get_voxels = function(latent_code,
                                                  resolution = 100,
                                                  sphere_only = TRUE) {

                              sequ <- seq(-1.1, 1.1, length.out = resolution)
                              ind <- seq_len(resolution)
                              indices <- expand.grid(ind, ind, ind)
                              pts <- expand.grid(sequ, sequ, sequ, KEEP.OUT.ATTRS = FALSE)

                              if(sphere_only) {
                                mask <- sqrt(rowSums(pts^2)) < 1.1
                                pts <- pts[mask, ]
                                indices <- indices[mask, ]
                              }

                              pts2 <- torch::torch_tensor(as.matrix(pts))
                              latent_codes <- latent_code$`repeat`(c(pts2$shape[1L], 1L))
                              sdf <- self$forward(pts2, latent_codes)

                              voxels <- array(1,
                                              dim = c(resolution,
                                                      resolution,
                                                      resolution))
                              voxels[as.matrix(indices)] <- as.array(sdf)

                              list(voxels = voxels, x = sequ, y = sequ, z = sequ)


                            },

                            get_mesh = function(latent_code, resolution = 100,
                                                smooth = FALSE) {

                              rlang::check_installed("rmarchingcubes")
                              rlang::check_installed("rgl")

                              voxels <-self$get_voxels(latent_code)
                              mesh <- rmarchingcubes::contour3d(voxels$voxels,
                                                                0,
                                                                x = voxels$x,
                                                                y = voxels$y,
                                                                z = voxels$z)
                              mesh2 <- rgl::mesh3d(vertices = rgl::asHomogeneous2(t(mesh$vertices)),
                                                   triangles = t(mesh$triangles),
                                                   normals = mesh$normals)

                              if(smooth) {
                                rlang::check_installed("Rvcg")
                                mesh2 <- Rvcg::vcgSmooth(mesh2)
                              }

                              mesh2

                            },

                            get_surface_points = function(latent_code, sample_size = 100000,
                                                          max_iter = 10,
                                                          sdf_cutoff = 0.0001,
                                                          ...) {

                              points <- get_points_in_unit_sphere(n = sample_size,
                                                                  device = self$parameters[[1]]$device)

                              latent_codes <- latent_code$`repeat`(c(points$shape[1L], 1L))
                              sdf <- self$forward(points, latent_codes)
                              keep_points <- points[sdf < sdf_cutoff, ]
                              points <- points[!sdf < sdf_cutoff, ]

                              iter <- 0
                              while(nrow(points) > 0 && iter < max_iter) {
                                iter <- iter + 1

                                points$requires_grad = TRUE
                                latent_codes <- latent_code$`repeat`(c(points$shape[1L], 1L))
                                sdf <- self$forward(points, latent_codes)
                                sdf$backward(torch::torch_ones_like(sdf, device = self$parameters[[1]]$device))

                                torch::with_no_grad({
                                  normals <- points$grad
                                  normals <- normals / torch::torch_norm(normals, dim = 2L)$unsqueeze(dim = 2L)
                                  points <- points$detach()
                                  points <- points - normals * sdf$unsqueeze(dim = 2L)
                                  sdf <- self$forward(points, latent_codes)
                                  keep_points <- torch::torch_cat(list(keep_points,
                                                                       points[sdf < sdf_cutoff, ]),
                                                                  dim = 1L)
                                  points <- points[!sdf < sdf_cutoff, ]
                                })
                                # Move points towards surface by the amount given by the signed distance

                              }


                              return(keep_points)

                            }


)

get_points_in_unit_sphere <- function(n, device = 'cpu') {

  x <- torch::torch_rand(as.integer(n * 2.5), 3L, device = device) * 2 - 1
  mask <- (torch::torch_norm(x, dim = 2L) < 1)$nonzero()$squeeze()
  mask <- mask[1:n]
  x <- x[mask, ]

  return(x)
}

map_cat <- function(.x, .f, ...) {
  purrr::map(.x, .f, ...) %>%
    torch::torch_cat()
}

