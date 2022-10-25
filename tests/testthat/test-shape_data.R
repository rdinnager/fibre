test_that("shape_data_inla() returns the correct data shapes", {
  ## 2 outcomes, 2 predictors, 2 pfcs, two latents
  dat <- hardhat::mold(lnBMR + lnGS ~ lnMass + endo + bre(phlo) + bre(phlo^2, latent = 2), 
                       data = vert_bmr,
                       blueprint = fibre_formula_blueprint())
  c(phyfs, rate_dists, hypers, latents) %<-% purrr::transpose(dat$extras$model_info)
  shaped <- shape_data_inla(phyfs, dat$predictors,
                            dat$outcomes, latents)
  expect_true(check_inla_dims(shaped$y, shaped$dat, shaped$A,
                                 2, 2, 2, 2))
  
  ## 2 outcomes, 2 predictors, 1 pfcs, two latents
  dat <- hardhat::mold(lnBMR + lnGS ~ lnMass + endo + bre(phlo, latent = 2), 
                       data = vert_bmr,
                       blueprint = fibre_formula_blueprint())
  c(phyfs, rate_dists, hypers, latents) %<-% purrr::transpose(dat$extras$model_info)
  shaped <- shape_data_inla(phyfs, dat$predictors,
                            dat$outcomes, latents)
  expect_true(check_inla_dims(shaped$y, shaped$dat, shaped$A,
                                 ny = 2, nlatent = 2, npredictors = 2, npfc = 1))
  
  ## 1 outcomes, 3 predictors, 2 pfcs, no latents
  dat <- hardhat::mold(lnBMR ~ lnMass + lnGS + endo + bre(phlo) + bre(phlo^2), 
                       data = vert_bmr,
                       blueprint = fibre_formula_blueprint())
  c(phyfs, rate_dists, hypers, latents) %<-% purrr::transpose(dat$extras$model_info)
  shaped <- shape_data_inla(phyfs, dat$predictors,
                            dat$outcomes, latents)
  expect_true(check_inla_dims(shaped$y, shaped$dat, shaped$A,
                                 ny = 1, nlatent = 0, npredictors = 3, npfc = 2))
  
})
