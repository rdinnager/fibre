prior_brownian <- function(phy, data, family, rate_order = c("first", "second")) {
  rate_order <- match.arg(rate_order)
  tree2 <- ape::makeNodeLabel(tree)
  phy_vcv <- MCMCglmm::inverseA(tree2, nodes = "ALL", scale = FALSE)$A
  ord_char <- c(tree2$tip.label, tree2$node.label[-1])
  phy_vcv <- phy_vcv[ord_char, ord_char]
  id <- 1:nrow(phy_vcv)
  mod <- INLA::inla(data ~ f(id, model = "generic0", Cmatrix = phy_vcv),
                    data = list(data = c(data, rep(NA, ape::Nnode(tree2) - 1))),
                    family = family,
                    control.compute = list(waic = TRUE))
  if(rate_order == "second") {
    av_rate <- 1 / mod$summary.hyperpar["Precision for id", ]$mean
  } else {
    edge <- tree2$edge
    rownames(edge) <-
    mod$summary.random$id
  }

}

prior_basic <- function(phy, dat, A_mat) {
  dat_sd <- sd(dat[ , resp])
  l <- A_mat > 0
  e_var <- (sqrt(dat_sd / (sum(A_mat[l]) / length(phy$tip.label)) /
                          (sum(l) / length(phy$tip.label)))) * 3
  e_var
}
