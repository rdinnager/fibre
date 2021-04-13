#' Make to the root-to-tip matrix of a phylogeny
#'
#' @param phy
#' @param return_nodes
#'
#' @return
#' @export
#'
#' @examples
make_L <- function(phy,
                   return_nodes = c("tips", "internal", "both"),
                   return_type = c("matrix", "list")) {

  return_nodes <- match.arg(return_nodes)
  return_type <- match.arg(return_type)

  temp_phy <- phy
  temp_phy$tip.label <- as.character(seq_along(temp_phy$tip.label))
  temp_phy$node.label <- as.character(length(temp_phy$tip.label) +
                                        seq_len(temp_phy$Nnode))


  ig <- igraph::as.igraph(temp_phy, directed = TRUE)

  if(!is.null(temp_phy$edge.length)) {
    edges_subtending <- match(as.numeric(names(igraph::V(ig))), temp_phy$edge[ , 2])
    igraph::vertex_attr(ig, "brlen") <- temp_phy$edge.length[edges_subtending]
    igraph::vertex_attr(ig, "brlen")[is.na(igraph::vertex_attr(ig, "brlen"))] <- 0
  }

  ## find root
  degs <- igraph::degree(ig, mode = "in")
  root <- names(degs)[degs == 0]

  n_nodes <- length(igraph::V(ig))

  if(return_nodes == "tips" | return_nodes == "both") {
    tips <- as.character(1:length(phy$tip.label))
  }

  if(return_nodes == "internal" | return_nodes == "both") {
    internal <- as.character((length(phy$tip.label) + 1):n_nodes)
  }

  if(return_nodes == "tips") {
    nodes <- tips
  }

  if(return_nodes == "internal") {
    nodes <- internal
  }

  if(return_nodes == "both") {
    nodes <- c(tips, internal)
  }

  root_to_tip <- igraph::shortest_paths(ig, from = root, to = nodes, mode = "out")

  rtp_mat <- lapply(root_to_tip$vpath, function(x) tabulate(x, n_nodes))
  rtp_mat <- do.call(cbind, rtp_mat)

  rtp_mat <- t(rtp_mat * igraph::vertex_attr(ig, "brlen"))

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

  rtp_mat

}
