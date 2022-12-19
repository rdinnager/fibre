library(phyf)
library(fibre)
library(torch)
library(gifski)
data(bird_beak_codes)

sdfnet <- load_bird_beak_model()

latent <- bird_beak_codes %>%
  dplyr::select(dplyr::starts_with("latent_")) %>%
  dplyr::slice(1500) %>%
  unlist() %>%
  matrix(nrow = 1) %>%
  torch_tensor(requires_grad = FALSE)

sdfnet$cuda()
#debug(sdfnet$get_normals)

test <- sdfnet$render_image(latent_code = latent,
                            camera_position = get_camera_position(-3, 150, 200),
                            resolution = 600,
                            iterations = 1000,
                            ssaa = 1,
                            threshold = 0.0005,
                            cuda = TRUE,
                            max_ray_move = 0.02,
                            verbose = TRUE)

mesh <- sdfnet$get_mesh(latent, resolution = 400)
rgl::shade3d(mesh, col = "red")

for(i in 1:180) {
  im <- try(sdfnet$render_image(latent_code = latent,
                            camera_position = get_camera_position(-2, i * 2, 179),
                            iterations = 1000,
                            ssaa = 1,
                            cuda = TRUE))

  if(inherits(im, "try-error")) {
    message("Error on image ", i)
  } else {
    imager::save.image(im, file.path("testing/images", paste0("test_", stringr::str_pad(i, 3, pad = "0"), ".png")))
  }

}

gif <- gifski::gifski(list.files("testing/images", full.names = TRUE), delay = 0.01)

sd1 <- load_state_dict(system.file("models/sdf_net.pt", package = "fibre"))
sdfnet <- sdf_net()
sdfnet$load_state_dict(sd1)

latent <- bird_beak_codes %>%
  dplyr::select(dplyr::starts_with("latent_")) %>%
  dplyr::slice(1500) %>%
  unlist() %>%
  torch_tensor()

bird_beak_codes$label[1500]

# sdfnet$cuda()
# latent <- latent$cuda()

sdfnet$eval()

pps <- sdfnet$get_surface_points(latent, 10000)
rgl::points3d(as.matrix(pps))

mesh <- sdfnet$get_mesh(latent)

rgl::shade3d(mesh, col = "red")

data3js <- r3js::plot3js(
  x = voxels$x,
  y = voxels$y,
  z = voxels$z,
  type = "n"
)

# Add shape according to the calculated contours
data3js <- r3js::shape3js(
  data3js,
  vertices = mesh$vertices,
  faces = mesh$triangles,
  normals = mesh$normals,
  col = "red"
)

# View the plot
r3js::r3js(data3js)


library(cgalMeshes)
test <- cgalMeshes::AFSreconstruction(as.array(pps))
rgl::shade3d(test$getMesh(), col = "red")

