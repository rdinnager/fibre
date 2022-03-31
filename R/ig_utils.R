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
#' @param order Can be `"first_order"` or `"second_order"`, specifying whether to return
#' the root-to-tip matrix in a first order or a second order specification. See documentation
#' of `fibre()` for a description of the difference.
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
                          order = c("first_order", "second_order"),
                          sparse = TRUE) {

  return_nodes <- match.arg(return_nodes)
  return_type <- match.arg(return_type)

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

  if(return_nodes == "tips" | return_nodes == "both" | return_ages | order == "second_order") {
    tips <- as.character(1:length(phy$tip.label))
  }

  if(return_nodes == "internal" | return_nodes == "both" | return_ages | order == "second_order") {
    internal <- as.character((length(phy$tip.label) + 1):n_nodes)
  }

  if(return_nodes == "tips") {
    nodes <- tips
  }

  if(return_nodes == "internal") {
    nodes <- internal
  }

  if(return_nodes == "both" | return_ages | order == "second_order") {
    nodes <- c(tips, internal)
  }

  root_to_tip <- igraph::shortest_paths(ig, from = root, to = nodes, mode = "out", output = "vpath")

  if(order == "second_order" | return_ages) {
    brlen <- igraph::vertex_attr(ig, "brlen", as.character(nodes))
    names(brlen) <- as.character(nodes)
    lens <- sapply(root_to_tip$vpath, function(x) c(sum(brlen[names(x[-length(x)])]),
                                                    brlen[names(x[length(x)])]))
    colnames(lens) <- names(brlen)
  }

  rtp_mat <- build_root2tip(root_to_tip, n_nodes, sparse = sparse)

  if(order == "second_order") {
    len_mat <- t(rtp_mat) * apply(lens, 2, sum)
    rtp_mat <- len_mat - t(rtp_mat * lens[1, names(igraph::V(ig))])
  } else {
    rtp_mat <- t(rtp_mat * igraph::vertex_attr(ig, "brlen"))
  }

  new_names <- c(phy$tip.label, temp_phy$node.label)

  colnames(rtp_mat) <- new_names[as.numeric(names(igraph::V(ig)))]
  rtp_mat <- rtp_mat[ , c(temp_phy$node.label, phy$tip.label)]

  if(return_nodes == "tips") {
    rownames(rtp_mat) <- phy$tip.label
  }
  if(return_nodes == "internal") {
    rownames(rtp_mat) <- temp_phy$node.label
  }
  if(return_nodes == "both") {
    rownames(rtp_mat) <- new_names
  }

  rtp_mat <- rtp_mat[ , -1]

  if(return_nodes == "both" & return_type == "list") {
    rtp_mat <- list(tips = rtp_mat[1:length(temp_phy$tip.label), ],
                    internal = rtp_mat[(length(temp_phy$tip.label) + 1):n_nodes, ])
  }

  if(return_ages) {
    attr(rtp_mat, "ages") <- lens[1, ]
  }

  rtp_mat

}

build_root2tip <- function(ig_out, n_nodes, sparse = TRUE) {

  if(!sparse) {

    rtp_mat <- lapply(ig_out$vpath, function(x) tabulate(x, n_nodes))
    do.call(cbind, rtp_mat)

  } else {

    js <- rep(seq_along(ig_out$vpath),
              lengths(ig_out$vpath))

    ig_out <- unlist(ig_out$vpath)

    Matrix::sparseMatrix(
      j = js,
      i = fastmatch::fmatch(ig_out, seq_len(n_nodes), nomatch = 0),
      x = 1
    )
  }

}
