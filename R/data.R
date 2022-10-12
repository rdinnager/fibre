#' AVONET Bird Trait Data with Phylogeny
#'
#' The AVONET Bird Trait Database joined to a `pf` object (for `{phyf}`)
#'
#' @format ## `who`
#' A 'pf' data frame (subclasses `tibble`) with 13,338 rows and 38 columns:
#' \describe{
#'   \item{node_name}{Node labels including species name for the tip labels}
#'   \item{phyf}{The phylogenetic flow column which stores the phylogenetic information}
#'   \item{Species3, Family3, Order3}{Taxonomic names -- Names of the species, family and order, respectively}
#'   \item{Total.individuals}{Number of individuals used to measure the data}
#'   \item{Beak.Length_Culmen:Species.Status}{Various traits of the bird species, see Source section to get more detailed information}
#'   ...
#' }
#' @source <https://figshare.com/articles/dataset/AVONET_morphological_ecological_and_geographical_data_for_all_birds_Tobias_et_al_2021_Ecology_Letters_/16586228>
"avonet"
