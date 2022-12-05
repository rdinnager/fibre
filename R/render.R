projection_matrix <- function() {
  matrix(
  c(1.73205081, 0,           0,           0,
    0,          1.73205081,  0,           0,
    0,          0,          -1.02020202, -0.2020202,
    0,          0,          -1,           0         ),
  nrow = 4, ncol = 4)
}

get_rotation_matrix <- function(angle, axis='y') {

  rlang::check_installed("orientlib")
  #rlang::check_installed("rgl")
  angles <- c(z = 0, y = 0, x = 0)
  angles[axis] <- angle * (pi/180)
  rot <- t(orientlib::rotmatrix(orientlib::eulerzyx(angles))@x[ , , 1])

  mat <- diag(4)
  mat[1:3, 1:3] <- rot
  return (mat)

}

get_camera_transform <- function(camera_distance, rotation_y, rotation_x = 0, project = FALSE) {

  #rlang::check_installed("rgl")

  camera_transform <- diag(4)
  camera_transform[3, 4] <- -camera_distance
  camera_transform <- camera_transform %*% get_rotation_matrix(rotation_x, axis = 'x')
  camera_transform <- camera_transform %*% get_rotation_matrix(rotation_y, axis = 'y')

  # camera_transform <- rgl::translationMatrix(0, 0, -camera_distance)
  # camera_transform <- camera_transform %*% get_rotation_matrix(rotation_x, axis = 'x')
  # camera_transform <- camera_transform %*% get_rotation_matrix(rotation_y, axis = 'y')

  if (project) {
    camera_transform <- projection_matrix() %*% camera_transform
  }

  return(camera_transform)
}

get_camera_position <- function(camera_distance, rotation_y, rotation_x = 0, project = FALSE) {
  camera_transform <- get_camera_transform(camera_distance, rotation_y, rotation_x, project)
  camera_position <- (solve(camera_transform) %*% matrix(c(0, 0, 0, 1), ncol = 1))[1:3, ]
  camera_position
}

default_camera <- function() {
  camera_transform <- get_camera_transform(2.2, 147, 20)
  camera_position <- (solve(camera_transform) %*% matrix(c(0, 0, 0, 1), ncol = 1))[1:3, ]
  camera_position
}

default_light <- function() {
  light_matrix <- get_camera_transform(6, 164, 50)
  light_position <- (solve(light_matrix) %*% matrix(c(0, 0, 0, 1), ncol = 1))[1:3, ]
  light_position
}

norm_vec <- function(x) sqrt(sum(x^2))

get_shadows <- function(sdf_net, points, light_position, latent_code = NULL, threshold = 0.001, sdf_offset=0, radius=1.0, cuda=FALSE, batch_size = 50000) {
  ray_directions <- t(t(-points) + light_position)
  ray_directions <- ray_directions / as.vector(as.matrix(torch::torch_norm(torch::torch_tensor(ray_directions), dim = -1)))

  points <- points + ray_directions * 0.1

  indices <- seq_len(nrow(points))
  shadows <- integer(nrow(points))

  for (i in seq_len(200)) {
    test_points <- torch::torch_tensor(points[indices, , drop = FALSE], dtype = torch::torch_float32())
    # if(latent_code$shape[1] == 1) {
    #   latent_code <- latent_code$`repeat`(c(test_points$shape[1], 1L))
    # }
    sdf <- sdf_net$evaluate_in_batches(test_points, latent_code, return_cpu = TRUE, cuda = cuda, batch_size = batch_size) + sdf_offset
    #sdf <- torch.clamp_(sdf, -0.1, 0.1)
    sdf <- as.vector(as.matrix(sdf))
    points[indices, ] <- points[indices, , drop = FALSE] + ray_directions[indices, ] * sdf

    hits <- (abs(sdf) < threshold)
    shadows[indices[hits]] <- 1
    indices <- indices[!hits]

    misses <- as.vector(as.matrix(((torch::torch_norm(torch::torch_tensor(points[indices, , drop = FALSE], dtype = torch::torch_float32()),
                                                      dim = 2) > radius)))) == 1
    indices <- indices[!misses]

    if (length(indices) < 2) {
      break
    }
    #print(length(indices))
  }

  shadows[indices] <- 1
  return(shadows)

}

render_image <- function(sdf_net, latent_code = NULL, resolution = 800, camera_position, light_position, threshold = 0.0005, sdf_offset = 0, iterations = 1000, ssaa = 2, radius = 1.0, crop = FALSE, color = c(R = 246 / 255, G = 236 / 255, B = 133 / 255), vertical_cutoff = NULL, max_ray_move = 0.05, plot = TRUE, cuda = FALSE, batch_size = 50000) {

  rlang::check_installed("Morpho")
  rlang::check_installed("einsum")
  rlang::check_installed("imager")
  message("started")

  camera_forward <- camera_position / norm_vec(camera_position) * -1

  camera_distance <- norm_vec(camera_position)
  up <- c(0, 1, 0)
  camera_right <- Morpho::crossProduct(camera_forward, up)
  camera_right <- camera_right / norm_vec(camera_right)
  camera_up <- Morpho::crossProduct(camera_forward, camera_right)
  camera_up <- camera_up / norm_vec(camera_up)

  screenspace_points <- tidyr::expand_grid(
    x = seq(-1, 1, length.out =  resolution * ssaa),
    y = seq(-1, 1, length.out = resolution * ssaa),
  ) %>%
    as.matrix()
  screenspace_points <- screenspace_points[ , c(2:1)]

  points <- matrix(camera_position, nrow = 1)[rep(1, nrow(screenspace_points)), ]
  #points = points.astype(np.float32)

  focal_distance <- 1.0 / tan(asin(radius / camera_distance))
  ray_directions <- screenspace_points[ , 1] * matrix(camera_right, nrow = 1)[rep(1, nrow(screenspace_points)), ] +
    screenspace_points[ , 2] * matrix(camera_up, nrow = 1)[rep(1, nrow(screenspace_points)), ] +
    focal_distance * matrix(camera_forward, nrow = 1)[rep(1, nrow(screenspace_points)), ]
  #ray_directions = ray_directions.transpose()
  ray_directions <- ray_directions / apply(ray_directions, 1, norm_vec)
  #ray_directions /= np.linalg.norm(ray_directions, axis=1)[:, np.newaxis]

  b <- einsum::einsum('ij,ij->i', points, ray_directions) * 2
  #b = np.einsum('ij,ij->i', points, ray_directions) * 2
  c = ((as.vector(camera_position) %*% as.vector(camera_position)) - radius * radius)[1, 1]
  distance_to_sphere = (-b - sqrt((b^2) - 4 * c)) / 2
  indices <- which(is.finite(distance_to_sphere))
  #indices = np.argwhere(np.isfinite(distance_to_sphere)).reshape(-1)

  points[indices, ] <- points[indices, ] + ray_directions[indices, ] * as.vector(distance_to_sphere[indices])

  #points <- torch::torch_tensor(points, dtype = torch::torch_float32())
  #ray_directions_t <- torch::torch_tensor(ray_directions, dtype = torch::torch_float32())

  #indices <- torch::torch_tensor(indices, dtype = torch::torch_int64())
  #model_mask <- torch::torch_zeros(points$shape[1], dtype = torch::torch_uint8())
  model_mask <- integer(nrow(points))

  for(i in seq_len(iterations)) {
    test_points <- torch::torch_tensor(points[indices, , drop = FALSE], dtype = torch::torch_float32())
    # if(latent_code$shape[1] == 1) {
    #   latent_code <- latent_code$`repeat`(c(test_points$shape[1], 1L))
    # }

    sdf <- sdf_net$evaluate_in_batches(test_points, latent_code, return_cpu = TRUE, cuda = cuda, batch_size = batch_size) + sdf_offset
    sdf <- torch::torch_clamp(sdf, -max_ray_move, max_ray_move)
    sdf <- as.vector(as.matrix(sdf))
    points[indices, ] <- points[indices, , drop = FALSE] + ray_directions[indices, , drop = FALSE ] * sdf

    #hits = abs(sdf) < threshold
    hits = sdf >= 0 & sdf < threshold
    model_mask[indices[hits]] <- 1
    indices = indices[!hits]

    misses <- as.vector(as.matrix(((torch::torch_norm(torch::torch_tensor(points[indices, , drop = FALSE], dtype = torch::torch_float32()),
                                dim = 2) > radius)))) == 1
    indices = indices[!misses]

    if (length(indices) < 2) {
      break
    }
    #print(range(sdf))
  }
  message("done raymarching")
  model_mask[indices] <- 1

  # test <- cbind(screenspace_points, model_mask) %>% as.data.frame()
  # ggplot(test, aes(x, y)) + geom_raster(aes(fill = model_mask)) + coord_equal() + theme_minimal()

  if(!is.null(vertical_cutoff)) {
    model_mask[points[ , 2] > vertical_cutoff] <- 0
    model_mask[points[ , 2] < -vertical_cutoff] <- 0
  }

  # if(latent_code$shape[1] == 1) {
  #     latent_code <- latent_code$`repeat`(c(points$shape[1], 1L))
  # }
  message("starting normals")
  print(latent_code$requires_grad)
  t_points <- torch::torch_tensor(points[model_mask == 1, , drop = FALSE], dtype = torch::torch_float32())
  normal = sdf_net$get_normals_in_batches(t_points, latent_code, cuda = cuda, batch_size = batch_size) %>%
    as.matrix()
  message("done normals")

  model_mask <- model_mask == 1
  model_points <- points[model_mask, ]

  # if(latent_code$shape[1] == 1) {
  #   latent_code <- latent_code$`repeat`(c(model_points$shape[1], 1L))
  # }
  seen_by_light = 1.0 - get_shadows(sdf_net, model_points, light_position, latent_code, radius = radius, sdf_offset = sdf_offset, cuda = cuda, batch_size = batch_size)
  message("done shadows")

  light_direction <- t(t(-model_points) + light_position)
  light_direction <- light_direction / apply(light_direction, 1, norm_vec)

  diffuse <- einsum::einsum('ij,ij->i', light_direction, normal)
  diffuse[diffuse < 0] <- 0
  diffuse[diffuse > 1] <- 1
  diffuse <- diffuse * seen_by_light

  reflect <- light_direction - as.vector(einsum::einsum('ij,ij->i', light_direction, normal)) * normal * 2
  reflect <- reflect / apply(reflect, 1, norm_vec)
  specular <- einsum::einsum('ij,ij->i', reflect, ray_directions[model_mask, ])
  specular[specular < 0] <- 0
  specular[specular > 1] <- 1
  specular <- (specular^20) * seen_by_light
  rim_light <- -einsum::einsum('ij,ij->i', normal, ray_directions[model_mask, ])
  rim_light[rim_light < 0] <- 0
  rim_light[rim_light > 1] <- 1
  rim_light <- 1.0 - rim_light
  rim_light = (rim_light^4) * 0.3

  color <- matrix(color, nrow = 1)[rep(1, length(diffuse)), ] * (as.vector(diffuse) * 0.5 + 0.5)
  color <- color + (as.vector(specular) * 0.3 + as.vector(rim_light))

  color[color < 0] <- 0
  color[color > 1] <- 1

  # test <- cbind(screenspace_points, col = "#000000") %>% as.data.frame() %>%
  #   mutate(x = as.numeric(x), y = as.numeric(y))
  # test[model_mask, "col"] <- rgb(color[ , 1], color[ , 2], color[ , 3])
  # ggplot(test, aes(x, y)) + geom_raster(aes(fill = col)) + coord_equal() +
  #   scale_fill_identity() +
  #   theme_minimal()

  bg <- c(1, 1, 1)
  col_mat <- matrix(rep(bg, nrow(screenspace_points)), nrow = nrow(screenspace_points))
  col_mat[model_mask, ] <- color
  img_arr <- array(t(col_mat), dim = c(3, resolution * ssaa, resolution * ssaa))
  img_arr <- aperm(img_arr, c(2, 3, 1))
  img <- imager::as.cimg(img_arr)

  message("done image")

  if(plot) {
    plot(img)
  }

  # ground_points = ray_directions[:, 1] < 0
  # ground_points[model_mask] = 0
  # ground_points = np.argwhere(ground_points).reshape(-1)
  # ground_plane = np.min(model_points[:, 1]).item()
  # points[ground_points, :] -= ray_directions[ground_points, :] * ((points[ground_points, 1] - ground_plane) / ray_directions[ground_points, 1])[:, np.newaxis]
  # ground_points = ground_points[np.linalg.norm(points[ground_points, ::2], axis=1) < 3]
  #
  # ground_shadows = get_shadows(sdf_net, points[ground_points, :], light_position, latent_code, sdf_offset=sdf_offset)
  #
  # pixels = np.ones((points.shape[0], 3))
  # pixels[model_mask] = color
  # pixels[ground_points] -= ((1.0 - 0.65) * ground_shadows)[:, np.newaxis]
  # pixels = pixels.reshape((resolution * ssaa, resolution * ssaa, 3))
  #
  #
  # if crop:
  #   from util import crop_image
  # pixels = crop_image(pixels, background=1)
  #
  # image = Image.fromarray(np.uint8(pixels * 255) , 'RGB')

  if (ssaa != 1) {
   img <- imager::resize(img, resolution, resolution, interpolation_type = 6)
  }

  return(img)
}
