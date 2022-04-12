test_that("make_tip2root works as expected", {
  set.seed(224235)
  test_tree <- rcoal(100)
  nnode <- Nnode(test_tree) - 1L ## subtract root node which is not included in root2tip matrix
  ntip <- Ntip(test_tree)

  rtp_1 <- make_root2tip(test_tree,
                          return_nodes = "tips",
                          return_type = "matrix",
                          return_ages = FALSE,
                          order = "first",
                          sparse = TRUE)

  expect_s4_class(rtp_1, "dgCMatrix")
  expect_identical(dim(rtp_1), c(ntip, ntip + nnode))
  expect_equal(var(Matrix::rowSums(rtp_1)), 0)

  rtp_2 <- make_root2tip(test_tree,
                          return_nodes = "internal",
                          return_type = "matrix",
                          return_ages = FALSE,
                          order = "first",
                          sparse = TRUE)

  expect_s4_class(rtp_2, "dgCMatrix")
  expect_identical(dim(rtp_2), c(nnode, ntip + nnode))
  expect_false(var(Matrix::rowSums(rtp_2)) == 0)

  expect_true(all(Matrix::rowSums(rtp_2) < mean(Matrix::rowSums(rtp_1))))

  rtp_3 <- make_root2tip(test_tree,
                          return_nodes = "both",
                          return_type = "matrix",
                          return_ages = FALSE,
                          order = "first",
                          sparse = TRUE)

  expect_s4_class(rtp_3, "dgCMatrix")
  expect_identical(dim(rtp_3), c(ntip + nnode, ntip + nnode))
  expect_equal(var(Matrix::rowSums(rtp_3[1:ntip, ])), 0)

  rtp_4 <- make_root2tip(test_tree,
                          return_nodes = "both",
                          return_type = "list",
                          return_ages = FALSE,
                          order = "first",
                          sparse = TRUE)

  expect_type(rtp_4, "list")
  expect_s4_class(rtp_4[[1]], "dgCMatrix")
  expect_s4_class(rtp_4[[2]], "dgCMatrix")
  expect_length(rtp_4, 2)
  expect_identical(dim(rtp_4[[1]]), c(ntip, ntip + nnode))
  expect_identical(dim(rtp_4[[2]]), c(nnode, ntip + nnode))
  expect_equal(rtp_1, rtp_4[[1]])
  expect_equal(rtp_2, rtp_4[[2]])
  expect_equal(var(Matrix::rowSums(rtp_4[[1]])), 0)

  rtp_5 <- make_root2tip(test_tree,
                          return_nodes = "both",
                          return_type = "list",
                          return_ages = FALSE,
                          order = "second",
                          sparse = TRUE)

  expect_type(rtp_5, "list")
  expect_s4_class(rtp_5[[1]], "dgCMatrix")
  expect_s4_class(rtp_5[[2]], "dgCMatrix")
  expect_length(rtp_5, 2)
  expect_identical(dim(rtp_5[[1]]), c(ntip, ntip + nnode))
  expect_identical(dim(rtp_5[[2]]), c(nnode, ntip + nnode))
  expect_identical(Matrix::nnzero(rtp_5[[1]]), Matrix::nnzero(rtp_1))
  expect_identical(Matrix::nnzero(rtp_5[[2]]), Matrix::nnzero(rtp_2))
  expect_failure(expect_equal(rtp_1, rtp_5[[1]]))
  expect_failure(expect_equal(rtp_2, rtp_5[[2]]))
  expect_equal(var(Matrix::rowSums(rtp_4[[1]])), 0)

  rtp_6 <- make_root2tip(test_tree,
                          return_nodes = "both",
                          return_type = "list",
                          return_ages = FALSE,
                          order = "both",
                          sparse = TRUE)

  expect_type(rtp_6, "list")
  expect_s4_class(rtp_6[[1]][[1]], "dgCMatrix")
  expect_s4_class(rtp_6[[1]][[2]], "dgCMatrix")
  expect_s4_class(rtp_6[[2]][[1]], "dgCMatrix")
  expect_s4_class(rtp_6[[2]][[2]], "dgCMatrix")
  expect_length(rtp_6, 2)
  expect_length(rtp_6[[1]], 2)
  expect_length(rtp_6[[2]], 2)
  expect_identical(dim(rtp_6[[1]][[1]]), c(ntip, ntip + nnode))
  expect_identical(dim(rtp_6[[1]][[2]]), c(nnode, ntip + nnode))
  expect_identical(dim(rtp_6[[2]][[1]]), c(ntip, ntip + nnode))
  expect_identical(dim(rtp_6[[2]][[2]]), c(nnode, ntip + nnode))
  expect_equal(rtp_1, rtp_6[[1]][[1]])
  expect_equal(rtp_2, rtp_6[[1]][[2]])


})
