## create AVONET dataset
library(phyf)
library(dplyr)

avonet_tree <- ape::read.nexus("extdata/HackettStage1_0001_1000_MCCTreeTargetHeights.nex")
avonet_dat <- readr::read_csv("extdata/AVONET3_BirdTree.csv")

avonet <- as_pf(avonet_tree)
avonet <- avonet %>%
  left_join(avonet_dat %>%
              mutate(node_name = gsub(" ", "_", Species3)))

usethis::use_data(avonet, overwrite = TRUE)
