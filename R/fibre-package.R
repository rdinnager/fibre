#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

.fibre_env <- new.env()
.fibre_env$multiple_order <- FALSE

.onAttach <- function(libname, pkgname) {

  if(!requireNamespace("INLA", quietly = TRUE)) {
    packageStartupMessage("Welcome to the fibre package. We noticed you don't have the R package INLA installed. INLA is required for fibre to work, so we recommend going to https://www.r-inla.org/download-install and following the instructions (since INLA is not on CRAN)")
  }

}

.onLoad <- function(libname, pkgname) {
  make_root2tip <<- memoise::memoise(make_root2tip)
  make_paths <<- memoise::memoise(make_paths)
  make_Cmatrix <<- memoise::memoise(make_Cmatrix)
  make_extraconstr <<- memoise::memoise(make_extraconstr)
}
