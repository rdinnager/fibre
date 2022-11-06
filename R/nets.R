#' A signed distance field based neural network model for generating 3d shapes
#'
#' @param n_latent Number of dimensions for the latent space
#' @param breadth Breadth of the multilayer perceptron networks
#'
#' @return A `torch::nn_module()`
#' @export
#'
#' @examples
#' sdf_net()
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
                              input <- torch::torch_cat(list(points, latent_codes), dim = 2L)
                              x <- self$layers1(input)
                              x <- torch::torch_cat(list(x, input), dim = 2L)
                              x <- self$layers2(x)
                              return(x$squeeze())
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

# class SDFNet(SavableModule):
#
#     def
#
#     def evaluate_in_batches(self, points, latent_code, batch_size=100000, return_cpu_tensor=True):
#         latent_codes = latent_code.repeat(batch_size, 1)
#         with torch.no_grad():
#             batch_count = points.shape[0] // batch_size
#             if return_cpu_tensor:
#                 result = torch.zeros((points.shape[0]))
#             else:
#                 result = torch.zeros((points.shape[0]), device=points.device)
#             for i in range(batch_count):
#                 result[batch_size * i:batch_size * (i+1)] = self.forward(points[batch_size * i:batch_size * (i+1), :], latent_codes)
#             remainder = points.shape[0] - batch_size * batch_count
#             result[batch_size * batch_count:] = self.forward(points[batch_size * batch_count:, :], latent_codes[:remainder, :])
#         return result
#
#     def get_voxels(self, latent_code, voxel_resolution, sphere_only=True, pad=True):
#         if not (voxel_resolution, sphere_only) in sdf_voxelization_helper:
#             helper_data = SDFVoxelizationHelperData(self.device, voxel_resolution, sphere_only)
#             sdf_voxelization_helper[(voxel_resolution, sphere_only)] = helper_data
#         else:
#             helper_data = sdf_voxelization_helper[(voxel_resolution, sphere_only)]
#
#         with torch.no_grad():
#             distances = self.evaluate_in_batches(helper_data.sample_points, latent_code).numpy()
#
#         if sphere_only:
#             voxels = np.ones((voxel_resolution, voxel_resolution, voxel_resolution), dtype=np.float32)
#             voxels[helper_data.unit_sphere_mask] = distances
#         else:
#             voxels = distances.reshape(voxel_resolution, voxel_resolution, voxel_resolution)
#             if pad:
#                 voxels = np.pad(voxels, 1, mode='constant', constant_values=1)
#
#         return voxels
#
#     def get_mesh(self, latent_code, voxel_resolution = 64, sphere_only = True, raise_on_empty=False, level=0):
#         size = 2 if sphere_only else 1.41421
#
#         voxels = self.get_voxels(latent_code, voxel_resolution=voxel_resolution, sphere_only=sphere_only)
#         voxels = np.pad(voxels, 1, mode='constant', constant_values=1)
#         try:
#             vertices, faces, normals, _ = skimage.measure.marching_cubes_lewiner(voxels, level=level, spacing=(size / voxel_resolution, size / voxel_resolution, size / voxel_resolution))
#         except ValueError as value_error:
#             if raise_on_empty:
#                 raise value_error
#             else:
#                 return None
#
#         vertices -= size / 2
#         mesh = trimesh.Trimesh(vertices=vertices, faces=faces, vertex_normals=normals)
#         return mesh
#
#     def get_normals(self, latent_code, points):
#         if latent_code.requires_grad or points.requires_grad:
#             raise Exception('get_normals may only be called with tensors that don\'t require grad.')
#
#         points.requires_grad = True
#         latent_codes = latent_code.repeat(points.shape[0], 1)
#         sdf = self.forward(points, latent_codes)
#         sdf.backward(torch.ones(sdf.shape[0], device=self.device))
#         normals = points.grad
#         normals /= torch.norm(normals, dim=1).unsqueeze(dim=1)
#         return normals
#
#     def get_surface_points(self, latent_code, sample_size=100000, sdf_cutoff=0.1, return_normals=False):
#         points = get_points_in_unit_sphere(n=sample_size, device=self.device)
#         points.requires_grad = True
#         latent_codes = latent_code.repeat(points.shape[0], 1)
#
#         sdf = self.forward(points, latent_codes)
#
#         sdf.backward(torch.ones((sdf.shape[0]), device=self.device))
#         normals = points.grad
#         normals /= torch.norm(normals, dim=1).unsqueeze(dim=1)
#         points.requires_grad = False
#
#         # Move points towards surface by the amount given by the signed distance
#         points -= normals * sdf.unsqueeze(dim=1)
#
#         # Discard points with truncated SDF values
#         mask = torch.abs(sdf) < sdf_cutoff
#         points = points[mask, :]
#         normals = normals[mask, :]
#
#         if return_normals:
#             return points, normals
#         else:
#             return points
#
#     def get_surface_points_in_batches(self, latent_code, amount = 1000):
#         result = torch.zeros((amount, 3), device=self.device)
#         position = 0
#         iteration_limit = 20
#         while position < amount and iteration_limit > 0:
#             points = self.get_surface_points(latent_code, sample_size=amount * 6)
#             amount_used = min(amount - position, points.shape[0])
#             result[position:position+amount_used, :] = points[:amount_used, :]
#             position += amount_used
#             iteration_limit -= 1
#         return result
#
#     def get_inception_score(self, sample_size=1000, latent_variance=1):
#         import inception_score
#         if not inception_score.available_for_points:
#             return 0
#         POINTCLOUD_SIZE = 1000
#         points = torch.zeros((sample_size * POINTCLOUD_SIZE, 3), device=self.device)
#         distribution = torch.distributions.normal.Normal(0, latent_variance)
#         latent_codes = distribution.sample([sample_size, LATENT_CODE_SIZE]).to(self.device)
#         for i in range(sample_size):
#             points[i * POINTCLOUD_SIZE:(i+1)*POINTCLOUD_SIZE, :] = self.get_surface_points_in_batches(latent_codes[i, :], amount=POINTCLOUD_SIZE)
#         return inception_score.inception_score_points(points, sample_size)

get_points_in_unit_sphere <- function(n, device = 'cpu') {

  x <- torch::torch_rand(as.integer(n * 2.5), 3L, device = device) * 2 - 1
  mask <- (torch::torch_norm(x, dim = 2L) < 1)$nonzero()$squeeze()
  mask <- mask[1:n]
  x <- x[mask, ]

  return(x)
}

