tree_split_into_epochs <- function(phy, interval = NULL, k = NULL) {

  
  nodelens <- Matrix::rowSums(root2node(phy, "all"))
  
  max_path <- max(nodelens)
  
  if(is.null(interval) && !is.null(k)) {
    intervals <- seq(0, max_path, length.out = k)
  } else {
    intervals <- seq(0, max_path, by = interval)
  }
  
  epochs <- ivs::iv(start = intervals[-length(intervals)], end = intervals[-1])
  edges <- ivs::iv(start = nodelens[phy$edge[ , 1]], end = nodelens[phy$edge[ , 2]])
  
  overlaps <- ivs::iv_locate_overlaps(edges, epochs, no_match = "drop") %>%
    dplyr::bind_cols(ivs::iv_align(edges, epochs, locations = .))
  
  colnames(overlaps) <- c("edge_num", "epoch_num", "edge_iv", "epoch_iv")
  
  overlaps <- overlaps %>%
    dplyr::mutate(intersect = ivs::iv_pairwise_intersect(edge_iv, epoch_iv),
                  edge_start = ivs::iv_start(edge_iv),
                  intersect_start = ivs::iv_start(intersect) - edge_start,
                  intersect_end = ivs::iv_end(intersect) - edge_start,
                  intersect_iv = ivs::iv(intersect_start, intersect_end)) 
  
  overlaps
  
  
}

