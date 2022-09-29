root2node <- function(phy, return_nodes = c("tips", "all"),
                      order = c("first", "second")) {

  return_nodes <- match.arg(return_nodes)
  order <- match.arg(order)

  if(return_nodes == "all") {
    phy <- add_internal_tips(phy)
  }
  
  rtp_bin <- root2tip_binary(phy)

  edge_ord <- rev(ape::postorder(phy))
  node_ord <- c(ape::Ntip(phy) + 1, phy$edge[edge_ord, 2])

  rtp_bin <- rtp_bin[node_ord, ]
  rtp_bin <- rtp_bin[-1, ]

  if(order == "second") {
    node_to_tip <- apply_sparse_fun_mult(rtp_bin,
                                         function(x) rev(cumsum(x)))
  }

  rtp_bin <- Matrix::t(rtp_bin)

  lens <- phy$edge.length[edge_ord]

  rtp <- rtp_bin %*% Matrix::Diagonal(length(lens), lens)

  if(order == "second") {
    rtp <- rtp * Matrix::t(node_to_tip)
  }

  rtp


}

root2tip_binary <- function(phy) {

  np <- ape::nodepath(phy)
  rtp <- build_rtp(np, nrow(phy$edge) + 1)

  rtp

}

build_rtp <- function(paths, n_nodes, sparse = TRUE) {

  if(!sparse) {

    rtp_mat <- lapply(paths, function(x) tabulate(x, n_nodes))
    do.call(rbind, rtp_mat)

  } else {

    js <- rep(seq_along(paths),
              lengths(paths))

    nj <- length(paths)

    ig_out <- unlist(paths)

    nn <- seq_len(n_nodes)

    Matrix::sparseMatrix(
      j = js,
      i = fastmatch::fmatch(ig_out, nn, nomatch = 0),
      x = 1,
      dims = c(n_nodes, nj),
    )
  }

}

apply_sparse_fun <- function(m, f, ...) {
  vapply(listCols(m), f, FUN.VALUE=0.0, ...)
}

apply_sparse_fun_mult <- function(m, f, ...) {

  new_x <- lapply(listCols(m), f, ...)
  m2 <- m
  m2@x <- unlist(new_x)
  m2

}

listCols <- function(m) {

  res <- split(m@x, findInterval(seq_len(Matrix::nnzero(m)), m@p, left.open = TRUE))
  res

}

split_i <- function(m) {
  nzc <- which(diff(m@p) != 0)
  tm <- list()
  tm[nzc] <- split(m@i, findInterval(seq_len(Matrix::nnzero(m)), m@p, left.open = TRUE))
  attr(tm, "nzc") <- nzc
  tm
}

find_terminals <- function(m) {
  tm <- split_i(m)
  if(length(attr(tm, "nzc")) > 0) {
    res <- integer(length(tm))
    res[attr(tm, "nzc")] <- vapply(tm[attr(tm, "nzc")], function(x) x[length(x)] + 1L, FUN.VALUE = 0L)
  } else {
    res <- vapply(tm, function(x) x[length(x)] + 1L, FUN.VALUE = 0L)
  }
  res
}

remove_terminal <- function(m) {

  drop_last <- function(x) {
    x[length(x)] <- 0
    x
  }
  new_m <- Matrix::drop0(apply_sparse_fun_mult(m, drop_last))
  new_m

}

get_internal_paths <- function(m) {

  node_ms <- list()
  nodes_done <- NULL
  nn <- Matrix::nnzero(m)
  new_m <- m

  i <- 1

  while(nn > 0) {
    new_m <- remove_terminal(new_m)
    new_tm <- find_terminals(new_m)
    keep_paths <- which(!duplicated(new_tm) & !new_tm %in% nodes_done)
    new_nodes <- new_tm[keep_paths]
    nodes_done <- union(nodes_done, new_nodes)
    node_ms[[i]] <- new_m[ , keep_paths]
    nn <- Matrix::nnzero(new_m)
    i <- i + 1
    #print(nn)
  }

  node_m <- do.call(cbind, node_ms)
  tm <- find_terminals(node_m)

  node_m <- node_m[ , order(tm)][ , -1]
  node_m


}

add_internal_tips <- function(phy) {
  node_labels <- length(phy$tip.label) + seq_len(phy$Nnode)
  new_phy <- phangorn::add.tips(phy, as.character(node_labels),
                                node_labels, 0)
}