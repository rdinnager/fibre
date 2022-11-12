#' Function to simulate continuous trait value histories on a phylogeny.
#'
#' @param phy A phylogenetic tree (phylo object) on which to simulate traits
#' @param rate_model The type of rate model for how rates of evolution evolve
#' on the phylogeny: "continuous" for continuous Brownian motion evolution of rates,
#' or "discrete" for evolution of rate "classes" across the phylogeny, using an mk model.
#' @param temp_trend_rates What temporal trend in rates should there be? A positive
#' number for an increase, and negative number for a decrease with the magnitude controlling
#' the strength of linear change. This trend is added to rates simulated under the rate_model.
#' @param rate_change If \code{rate_model} is "continuous", this should be a single
#' positive number controlling how fast rates change continuously along the tree. If
#' \code{rate_model} is "discrete", this should be a transition matrix for the rate classes.
#' Or, if \code{rate_model} is "discrete", and this can be a length 2 numeric vector specifying
#' @param rates Only used if \code{rate_model} is "discrete", in which case this should
#' be a named vector whose values are the rates in each rate class, and whose names are the
#' rate class states (e.g. \code{c("1" = 3, "2" = 10)}). See \code{\link[phylotools]{sim.history}},
#' for more detail on how the discrete model works. Or, if an unnamed numeric vector of length
#' two, a mean and standard deviation parameterizing a normal distribution from which to
#' draw rates for each rate class. If \code{NULL}, rates will be drawn from a normals distribution
#' with mean = 0 and sd = 1.
#' @param anc Value of the trait at the root ancestor. For \code{rate_model = "discrete"},
#' can be a length 1 named vector where the name is the ancestral state, and the value is
#' the trait starting value. For \code{rate_model = "continuous"}, any names are ignored,
#' but should be length 2, where the first element is the ancestral trait value and the
#' second element is the ancestral rate of evolution.
#' @param internal Logical value. If \code{TRUE} return trait states at internal nodes.
#' @param nsim Number of simulation to run.
#' @param pos_strat ?
#' @param temp_trend_mean A temporal trend in rates.
#'
#' @return A vector or matrix (for \code{nsim > 1}) containing simulated trait values
#' for each tip if \code{internal = FALSE}, or for each node if \code{internal = TRUE}
#' @export
simulate_traits <- function(phy, rate_model = c("continuous", "discrete"), temp_trend_rates = 0,
                            rate_change, rates = NULL, anc = c("1" = 0), internal = FALSE,
                            nsim = 1, pos_strat = c("none", "log", "add_const"),
                            temp_trend_mean = 0) {

  rate_model <- match.arg(rate_model)

  if(temp_trend_rates != 0 | temp_trend_mean != 0) {
    node_times <- ape::node.depth.edgelength(phy)
    if(temp_trend_rates != 0) {
      node_rate_add <- node_times * temp_trend_rates
      edge_rate_add <- apply(phy$edge, 1, function(y) mean(node_rate_add[y])) * phy$edge.length
    }
    if(temp_trend_mean != 0) {
      node_add <- node_times * temp_trend_mean
    }
  }

  if(rate_model == "discrete") {

    if(!is.numeric(rates)) {
      stop("rates argument must be numeric.")
    }

    if(is.null(rates)) {
      warning("No rates specified. Drawing random values from a standard Normal distribution.")
      rates <- setNames(rnorm(nrow(rate_change), 0, 1), 1:nrow(rate_change))
    } else {
      if(is.null(names(rates)) & length(rates) == 2) {
        rates <- setNames(rnorm(nrow(rate_change), rates[1], rates[2]), 1:nrow(rate_change))
      }
    }

    ## simulate discrete character history
    tree <- sim.history(phy, rate_change, anc = names(anc)[1], nsim = nsim)

  } else {
    if(length(anc) < 2) {
      warning("For continuous model anc should be length 2: Assuming ancestral rate
              value is zero.")
      anc <- c(anc, 0)
    }
    if(length(anc) > 2) {
      warning("For continuous model anc should be length 2: Only using first 2 values.")
      anc <- anc[1:2]
    }
    rate_vals <- phytools::fastBM(phy, a = anc[2], sig2 = rate_change, internal = TRUE,
                                  nsim = nsim)

    if(nsim > 1) {
      edge_rates <- lapply(rate_vals,
                           function(z) apply(phy$edge, 1,
                                             function(y) mean(z[y])) * phy$edge.length)
      mapped.edge <- lapply(edge_rates,
                            function(y) matrix(y, ncol = 1,
                                               dimnames =
                                                 list(apply(phy$edge, 1, paste, collapse = ","),
                                                      "1")))
      tree <- lapply(mapped.edge, function(y) {tt <- phy; tt$mapped.edge <- y; tt})
    } else {

      tree <- phy

      edge_rates <- apply(phy$edge, 1, function(y) mean(rate_vals[y])) * phy$edge.length
      tree$mapped.edge <- matrix(edge_rates, ncol = 1,
                                 dimnames = list(apply(phy$edge, 1, paste, collapse = ","),
                                                 "1"))
    }

    rates <- c("1" = 1)

  }

  if(temp_trend_rates != 0) {
    if(nsim > 1) {
      tree <- lapply(tree, function(y) {
        y$mapped.edge <- cbind(y$mapped.edge,
                               edge_rate_add)
        nstate <- ncol(y$mapped.edge)
        colnames(y$mapped.edge)[nstate] <- nstate
        y
      })
      nstate <- ncol(tree[[1]]$mapped.edge)
    } else {
      tree$mapped.edge <- cbind(tree$mapped.edge,
                                edge_rate_add)
      nstate <- ncol(tree$mapped.edge)
      colnames(tree$mapped.edge)[nstate] <- nstate
    }
    rates <- c(rates, 1)
    names(rates)[nstate] <- nstate
  }

  if(pos_strat == "log") {
    tree$mapped.edge <- exp(tree$mapped.edge)
  }
  if(pos_strat == "add_const") {
    if(any(tree$mapped.edge < 0)) {
      tree$mapped.edge <- tree$mapped.edge + abs(min(tree$mapped.edge)) + 0.0001
    }
  }

  ## simulate traits
  if(nsim > 1) {
    x <- lapply(tree, function(y) sim.rates(y, rates, internal = internal, anc = anc[1]))
    x <- do.call(cbind, x)
    colnames(x) <- paste0("sim_", seq_len(ncol(x)))
  } else {
    x <- sim.rates(tree, rates, internal = internal, anc = anc[1])
  }

  if(temp_trend_mean != 0) {
    x <- x + node_add
  }

  x
}
