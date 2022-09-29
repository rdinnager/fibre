vcv_fast <- function(phy, internal = FALSE) {
  root2node <- root2node(phy)
  vcv <- sqrt(Matrix::tcrossprod(root2node))
}