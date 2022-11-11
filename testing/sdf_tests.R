library(phyf)
library(fibre)
library(torch)
data(bird_beak_codes)

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

