#' Make to the root-to-tip matrix of a phylogeny
#'
#' @param phy An object of class `phylo` containing the phylogeny to generate a root-to-tip
#' matrix for. Must be a rooted phylogeny but can contain polytomies (that is `is.binary(phy)` can
#' be `FALSE`).
#' @param return_nodes Which nodes to return rows for? (Can be `"tips"`, `"internal"` or `"both"`).
#' @param return_type If `return_nodes == "both"`, pass `"matrix"` to return a `matrix` with
#' both node types bound together using `rbind`, or pass `"list"` to return two separate matrices
#' in a list (with `"tips"` as the first element and `"internal"` as the second). Ignored if
#' `return_nodes != "both"`.
#' @param return_ages Should a vector of node 'ages' be added as an attribute named `"ages"` to
#' the return value? 'ages' is defined as the total distance to the node from the phylogenies
#' root. Note this is only really an age if the phylogeny is a time tree and that it is the
#' opposite of the usual meaning of age, that is, it doesn't measure the age of the clade
#' descending from the node, but how much time (or change accumulated) since the common ancestor
#' of the whole phylogeny (the root) when this node is reached. This definition of 'age'
#' makes sense whether the phylogeny is a time-tree or not. If a time-tree the more usual
#' definition of age can be calculated by subtracted the focal node's age from the age of
#' the tips (which should all be the same for a time-tree).
#' @param return_parent If `TRUE`, the parent node, matching the root to tip matrix columns,
#' will be returned as an attribute `"parents"`.
#' @param order Can be `"first_order"` or `"second_order"`, specifying whether to return
#' the root-to-tip matrix in a first order or a second order specification. See documentation
#' of [fibre()] for a description of the difference.
#' @param sparse Should the returned root-to-tip matrix use a sparse matrix representation?
#' If `FALSE` a regular dense matrix will be returned. The sparse matrix representation is
#' recommended, especially for large trees. Beyond a certain size of tree, only a
#' sparse return will work, a dense matrix will exhaust all memory. Most root-to-tip
#' matrices are quite sparse, and so can represent very large trees where dense
#' matrices cannot (how sparse exactly depends on the shape of the tree).
#'
#' @return A `matrix` if `sparse = FALSE`, or a `Matrix::dgCMatrix` otherwise.
#' @export
#'
#' @importFrom Matrix t
#'
#' @examples
#' test_tree <- ape::rtree(1000)
#' rtp_mat <- make_root2tip(test_tree)
#' rtp_mat
make_root2tip <- function(phy,
                          return_nodes = c("tips", "internal", "both"),
                          return_type = c("matrix", "list"),
                          return_ages = FALSE,
                          return_parents = TRUE,
                          order = "first",
                          sparse = TRUE,
                          threads = 1) {

  return_nodes <- match.arg(return_nodes)
  return_type <- match.arg(return_type)
  order <- match.arg(order, c("first", "second", "both"),
                     several.ok = TRUE)

  temp_phy <- phy
  temp_phy$tip.label <- as.character(seq_along(temp_phy$tip.label))
  temp_phy$node.label <- as.character(length(temp_phy$tip.label) +
                                        seq_len(temp_phy$Nnode))


  ig <- igraph::as.igraph(temp_phy, directed = TRUE)

  if(!is.null(temp_phy$edge.length)) {
    edges_subtending <- fastmatch::fmatch(as.numeric(names(igraph::V(ig))), temp_phy$edge[ , 2])
    igraph::vertex_attr(ig, "brlen") <- temp_phy$edge.length[edges_subtending]
    igraph::vertex_attr(ig, "brlen")[is.na(igraph::vertex_attr(ig, "brlen"))] <- 0
  }

  ## find root
  degs <- igraph::degree(ig, mode = "in")
  root <- names(degs)[degs == 0]

  n_nodes <- length(igraph::V(ig))

  tips <- as.character(1:length(phy$tip.label))
  internal <- as.character((length(phy$tip.label) + 1):n_nodes)

  # if(return_nodes == "tips" || return_nodes == "both" || return_ages || order == "second" || order == "both") {
  #   tips <- as.character(1:length(phy$tip.label))
  # }
  #
  # if(return_nodes == "internal" || return_nodes == "both" || return_ages || order == "second" || order == "both") {
  #   internal <- as.character((length(phy$tip.label) + 1):n_nodes)
  # }
  #
  # if(return_nodes == "tips") {
  #   nodes <- tips
  # }
  #
  # if(return_nodes == "internal") {
  #   nodes <- internal
  # }
  #
  # if(return_nodes == "both" || return_ages || order == "second" || order == "both") {
  #   nodes <- c(tips, internal)
  # }

  if(threads > 1) {
    to <- c(tips, internal)
    splitter <- gl(threads, ceiling(length(to) / threads), length = length(to))
    groups <- split(to, splitter)
    cl <- parallel::makeCluster(threads)
    parallel::clusterExport(cl, c("ig", "root"), envir = environment())
    root_to_tip <- parallel::parLapply(groups,
                                  function(ti) igraph::shortest_paths(ig, from = root,
                                                                      to = ti,
                                                                      mode = "out",
                                                                      output = "vpath"),
                                  cl = cl)
    parallel::stopCluster(cl)
    root_to_tip <- list(vpath = do.call(c, lapply(root_to_tip, function(x) x$vpath)),
                        epath = NULL,
                        predecessors = NULL,
                        inbound_edges = NULL)

  } else {
    root_to_tip <- make_paths(ig, from = root, to = c(tips, internal))
  }

  if(order == "second" || order == "both" || return_ages) {
    brlen <- igraph::vertex_attr(ig, "brlen", as.character(c(tips, internal)))
    names(brlen) <- as.character(c(tips, internal))
    lens <- sapply(root_to_tip$vpath, function(x) c(sum(brlen[names(x[-length(x)])]),
                                                    brlen[names(x[length(x)])]))
    colnames(lens) <- names(brlen)
  }

  rtp_mat <- build_root2tip(root_to_tip, n_nodes, sparse = sparse)

  new_names <- c(phy$tip.label, temp_phy$node.label)

  # if(return_nodes == "tips") {
  #   rtp_mat_r <- temp_phy$tip.label
  # }
  # if(return_nodes == "internal") {
  #   rtp_mat_r <- temp_phy$node.label
  # }
  #if(return_nodes == "both") {
    rtp_mat_r <- new_names
  #}

  if(order == "second" || order == "both") {
    len_mat <- t(rtp_mat) * apply(lens, 2, sum)
    rtp_mat_2 <- len_mat - t(rtp_mat * lens[1, names(igraph::V(ig))])
    colnames(rtp_mat_2) <- new_names[as.numeric(names(igraph::V(ig)))]
    rtp_mat_2 <- rtp_mat_2[ , c(temp_phy$node.label, phy$tip.label)]
    rownames(rtp_mat_2) <- rtp_mat_r
    rtp_mat_2 <- rtp_mat_2[ , -1]
  }

  if(order == "first" || order == "both") {
    rtp_mat_1 <- t(rtp_mat * igraph::vertex_attr(ig, "brlen"))
    colnames(rtp_mat_1) <- new_names[as.numeric(names(igraph::V(ig)))]
    rtp_mat_1 <- rtp_mat_1[ , c(temp_phy$node.label, phy$tip.label)]
    rownames(rtp_mat_1) <- rtp_mat_r
    rtp_mat_1 <- rtp_mat_1[ , -1]
  }

  if(order == "first") {
    rtp_mat <- rtp_mat_1
  }
  if(order == "second") {
    rtp_mat <- rtp_mat_2
  }

  # if(return_nodes == "internal" || return_nodes == "both") {
  #   if(order == "both") {
  #     rtp_mat_1 <- rtp_mat_1[rownames(rtp_mat_1) != root, ]
  #     rtp_mat_2 <- rtp_mat_2[rownames(rtp_mat_2) != root, ]
  #   } else {
  #     rtp_mat <- rtp_mat[rownames(rtp_mat) != root, ]
  #   }
  # }

  if(return_nodes == "both" && return_type == "list") {
    if(order != "both") {
      rtp_mat <- list(tips = rtp_mat[1:length(temp_phy$tip.label), ],
                      internal = rtp_mat[(length(temp_phy$tip.label) + 1):(n_nodes), ])
    } else {
      rtp_mat <- list(first = list(tips = rtp_mat_1[1:length(temp_phy$tip.label), ],
                                   internal = rtp_mat_1[(length(temp_phy$tip.label) + 1):(n_nodes), ]),
                      second = list(tips = rtp_mat_2[1:length(temp_phy$tip.label), ],
                                    internal = rtp_mat_2[(length(temp_phy$tip.label) + 1):(n_nodes), ]))
    }
  } else {
    if(order == "both") {
      if(return_nodes == "tips") {
        rtp_mat_1 <- rtp_mat_1[1:length(temp_phy$tip.label), ]
        rtp_mat_2 <- rtp_mat_2[1:length(temp_phy$tip.label), ]
      }
      if(return_nodes == "internal") {
        rtp_mat_1 <- rtp_mat_1[(length(temp_phy$tip.label) + 1):(n_nodes), ]
        rtp_mat_2 <- rtp_mat_2[(length(temp_phy$tip.label) + 1):(n_nodes), ]
      }
      rtp_mat = list(first = rtp_mat_1,
                     second = rtp_mat_2)
    } else {
      if(return_nodes == "tips") {
        rtp_mat <- rtp_mat[1:length(temp_phy$tip.label), ]
      }
      if(return_nodes == "internal") {
        rtp_mat <- rtp_mat[(length(temp_phy$tip.label) + 1):(n_nodes), ]
      }
    }
  }

  if(return_ages) {
    ages <- lens[1, ]
    ages <- ages[names(ages) != root]
    attr(rtp_mat, "ages") <- ages
  }

  if(return_parents) {
    parents <- matrix(fastmatch::fmatch(temp_phy$edge,
                                        c(as.numeric(temp_phy$node.label),
                                          as.numeric(temp_phy$tip.label))),
                      nrow = nrow(temp_phy$edge),
                      ncol = ncol(temp_phy$edge))
    attr(rtp_mat, "parents") <- parents
  }

  cols2edges <- match_cols2edges(colnames(rtp_mat), phy)
  attr(rtp_mat, "cols2edges") <- cols2edges

  #class(rtp_mat) <- c(class(rtp_mat), "fibre_root2tip")

  rtp_mat

}

build_root2tip <- function(ig_out, n_nodes, sparse = TRUE) {

  if(!sparse) {

    rtp_mat <- lapply(ig_out$vpath, function(x) tabulate(x, n_nodes))
    do.call(cbind, rtp_mat)

  } else {

    js <- rep(seq_along(ig_out$vpath),
              lengths(ig_out$vpath))

    nj <- length(ig_out$vpath)

    ig_out <- unlist(ig_out$vpath)

    Matrix::sparseMatrix(
      j = js,
      i = fastmatch::fmatch(ig_out, seq_len(n_nodes), nomatch = 0),
      x = 1,
      dims = c(n_nodes, nj),
    )
  }

}

make_paths <- function(ig, from, to) {
  igraph::shortest_paths(ig, from = from, to = to, mode = "out", output = "vpath")
}

match_cols2edges <- function(col_names, phy) {
  tips <- which(col_names %in% phy$tip.label)
  col_names[tips] <- as.character(match(col_names[tips], phy$tip.label))
  nodes <- as.integer(col_names)
  edges <- match(nodes, phy$edge[ , 2])
  edges
}

find_nearest_node_dist <- function(phy, from, to, use = c("steps", "length", "both")) {

  use <- match.arg(use)

  temp_phy <- phy
  temp_phy$tip.label <- as.character(seq_along(temp_phy$tip.label))
  temp_phy$node.label <- as.character(length(temp_phy$tip.label) +
                                        seq_len(temp_phy$Nnode))


  ig <- igraph::as.igraph(temp_phy, directed = TRUE)

  if(!is.null(temp_phy$edge.length)) {
    edges_subtending <- fastmatch::fmatch(as.numeric(names(igraph::V(ig))), temp_phy$edge[ , 2])
    igraph::vertex_attr(ig, "brlen") <- temp_phy$edge.length[edges_subtending]
    igraph::vertex_attr(ig, "brlen")[is.na(igraph::vertex_attr(ig, "brlen"))] <- 0
  }

  ## find root
  # degs <- igraph::degree(ig, mode = "in")
  # root <- names(degs)[degs == 0]

  get_closest_dist <- function(tipward, rootward) {
    if(tipward <= rootward) {
      dist <- tipward
    } else {
      dist <- -rootward
    }
    dist
  }

  if(use == "steps" || use == "both") {
    tipward_dists <- igraph::distances(ig, as.character(from), as.character(to), mode = "out")
    rootward_dists <- igraph::distances(ig, as.character(from), as.character(to), mode = "in")
    tipward_closest <- apply(tipward_dists, 1, which.min)
    rootward_closest <- apply(rootward_dists, 1, which.min)
    tipward_closest <- colnames(tipward_dists)[tipward_closest]
    rootward_closest <- colnames(rootward_dists)[rootward_closest]
    tipward_step_dist <- apply(tipward_dists, 1, min)
    rootward_step_dist <- apply(rootward_dists, 1, min)

    closest_dist <- mapply(get_closest_dist, tipward_step_dist, rootward_step_dist)
  }

  if(use == "length" || use == "both") {

    weights <- igraph::vertex_attr(ig, "brlen")[igraph::tail_of(ig, igraph::E(ig))]

    tipward_dists_l <- igraph::distances(ig, as.character(from), as.character(to), mode = "out", weights = weights)
    rootward_dists_l <- igraph::distances(ig, as.character(from), as.character(to), mode = "in", weights = weights)
    tipward_closest_l <- apply(tipward_dists_l, 1, which.min)
    rootward_closest_l <- apply(rootward_dists_l, 1, which.min)
    tipward_closest_l <- colnames(tipward_dists_l)[tipward_closest_l]
    rootward_closest_l <- colnames(rootward_dists_l)[rootward_closest_l]
    tipward_len_dist <- apply(tipward_dists_l, 1, min)
    rootward_len_dist <- apply(rootward_dists_l, 1, min)

    closest_dist_l <- mapply(get_closest_dist, tipward_len_dist, rootward_len_dist)

  }

  if(use == "steps") {
    res <- closest_dist
  }

  if(use == "length") {
    res <- closest_dist_l
  }

  if(use == "both") {
    res <- list(steps = closest_dist,
                length = closest_dist_l)
  }

  res

}
