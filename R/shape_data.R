shape_data_inla <- function(pfcs, predictors,
                            outcomes, latents) {
  
  # return(list(pfcs = pfcs, predictors = predictors,
  #             outcomes = outcomes, latents = latents))
  
  ny <- ncol(outcomes)
  ylen <- nrow(outcomes)
  
  latent_df <- purrr::map(seq_along(pfcs), #pfcs, latents,
                            ~ make_latent_data(pfcs[[.x]], latents[[.x]], ylen, .x)) 
  latent_df <- expand_empty(!!!latent_df) %>%
    dplyr::bind_cols()
  
  latent_copy_df <- purrr::map(seq_along(pfcs), #latents,
                            ~ make_latent_copy(pfcs[[.x]], latents[[.x]], ny, ylen, .x))
  latent_copy_df <- expand_empty(!!!latent_copy_df) %>%
    dplyr::bind_cols()

  y_df <- purrr::map(seq_along(pfcs), #latents,
                      ~ make_y_data(pfcs[[.x]], latents[[.x]], ny, ylen, .x))
  y_df <-  expand_empty(!!!y_df) %>%
    dplyr::bind_cols()
  
  if("(Intercept)" %in% colnames(predictors)) {
    predictors <- dplyr::bind_cols(predictors %>%
                                     dplyr::select(-`(Intercept)`),
                                   predictors %>%
                                     dplyr::select(`(Intercept)`) %>%
                                     list() %>%
                                     rep( ny) %>%
                                     purrr::imap(function(x, y) {
                                       colnames(x) <- paste0(colnames(x), "_", y)
                                       x
                                     }) %>%
                                     dplyr::bind_cols())
  }
  
  dat_pred <- purrr::imap(predictors,
                          compress_data)
  
  dat_pred <- purrr::transpose(dat_pred)
  
  x <- dplyr::bind_cols(purrr::map(dat_pred$data,
                                   ~ .x[1, ])) %>%
    dplyr::slice(0)
  x <- dplyr::bind_rows(c(list(x), dat_pred$data))
  
  pred_A <- do.call(cbind, dat_pred$A)
  rownames(pred_A) <- paste0("y_", seq_len(nrow(pred_A)))
  colnames(pred_A) <- paste0("x_", seq_len(ncol(pred_A)))
  x_pfc <- phyf::pf_as_pfc(pred_A, is_tip = rep(TRUE, nrow(pred_A)))
  x_df <- dplyr::tibble(x_pfc = x_pfc)
  
  x_df <- dplyr::bind_rows(rep(list(x_df), ny))
  
  # latent_copy <- latent_df %>%
  #   dplyr::mutate(dplyr::across(.fns = ~ phyf::pf_ones(.x)))
  #colnames(latent_copy_df) <- paste0("copy_", colnames(latent_copy_df))
  
  x_all_y <- dplyr::bind_cols(expand_empty(x_df, y_df, latent_copy_df))
  x_all_A <- dplyr::bind_rows(x_all_y, latent_df)
  
  y_block <- as.data.frame(diag(seq_len(ny)))
  ys <- tibble_block(y_block, purrr::imap(outcomes, 
                                          ~ dplyr::tibble("{.y}" := .x)),
                     glue_names = FALSE)
  
  ys_latent <- purrr::map(latents,
                          ~ make_latent_y(.x, nrow(outcomes))) %>%
    dplyr::bind_rows()
  
  ys_all <- dplyr::bind_rows(ys, ys_latent)
  
  y_df_A <- purrr::map(y_df, phyf::pf_as_sparse)
  y_df_new <- purrr::imap_dfr(y_df_A,
                               ~ dplyr::tibble("{.y}_indexes" := seq_len(ncol(.x))))
  
  x_df_A <- purrr::map(x_df, phyf::pf_as_sparse)
  x_df_new <- purrr::imap_dfr(x_df_A,
                               ~ dplyr::tibble("{.y}_indexes" := seq_len(ncol(.x)))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(x_pfc_indexes = list(x[x_pfc_indexes, ])) %>%
    tidyr::unnest(cols = dplyr::everything())
  
  if(sum(unlist(latents)) > 0) {
    latent_copy_A <- purrr::map(latent_copy_df, phyf::pf_as_sparse)
    latent_copy_new <- purrr::imap_dfr(latent_copy_A, 
                                   ~ dplyr::tibble("{.y}_indexes" := seq_len(ncol(.x)),
                                                   "names" := paste0(.y, ":", colnames(.x)))) %>%
      tidyr::separate(.data$names, c("type", "latent", "ename"), sep = ":") %>%
      dplyr::select(-dplyr::all_of(c("type", "ename"))) %>%
      dplyr::group_by(.data$latent) %>%
      dplyr::group_nest() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(data = list(data %>%
                      stats::setNames(paste(latent, colnames(data), sep = ":")))) %>%
      dplyr::ungroup() %>%
      tidyr::unnest(dplyr::everything()) %>%
      dplyr::select(-dplyr::all_of("latent"))
    
    latent_df_A <- purrr::map(latent_df, phyf::pf_as_sparse)
    latent_df_new <- purrr::imap_dfr(latent_df_A,
                                     ~ dplyr::tibble("{.y}_indexes" := seq_len(ncol(.x))))
  } else {
    latent_copy_A <- MatrixExtra::emptySparse(format = "C")
    latent_copy_new <- dplyr::tibble()
    latent_df_A <- MatrixExtra::emptySparse(format = "C")
    latent_df_new <- dplyr::tibble()
  }
  
  #x_all_y <- dplyr::bind_cols(x_df, latent_copy)
  x_all <- dplyr::bind_rows(x_df_new, y_df_new, latent_copy_new, latent_df_new)
  
  A_all <- do.call(cbind, purrr::map(x_all_A, phyf::pf_as_sparse))

  check_inla_dims(ys_all, x_all, A_all, ny, sum(unlist(latents)), ncol(predictors), length(pfcs))
    
  return(list(dat = x_all, y = ys_all, A = A_all))
  
}

compress_data <- function(predictor, name) {

  compressed <- dplyr::tibble(predictor) %>%
    stats::setNames(name) %>%
    dplyr::group_by(.data[[name]]) %>%
    dplyr::group_data() %>%
    dplyr::mutate(num = seq_len(dplyr::n())) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(.cols = list(rep(num, length(.rows)))) %>%
    dplyr::ungroup()
  
  data <- compressed[1]
  A <- Matrix::sparseMatrix(unlist(compressed$.rows), 
                            unlist(compressed$.cols),
                            x = 1)
  list(data = data, A = A)
}

make_latent_data <- function(pfc, latent, ylen, num) {
  
  if(latent > 0) {
    d <- Matrix::.sparseDiagonal(latent, shape = "g")
    colnames(d) <- paste0("latent_", seq_len(latent)) 
    new_pfc <- phyf::pf_kronecker(d, pfc)
    # y_latent <- tibble_block(as.data.frame(as.matrix(d)),
    #                          list(dplyr::tibble(y = rep(NA, length(pfc)))))
    new_pfc <- dplyr::tibble("latent_pfc_{num}" := new_pfc)
    # y_latent <- dplyr::bind_cols(y_latent, dplyr::tibble(latent_pfc = new_pfc))
  } else {
    new_pfc <- dplyr::tibble()
  }
  
  new_pfc
  
}

make_latent_copy <- function(pfc, latent, ny, ylen, num) {
  
  if(latent > 0) {
    d <- matrix(1, ncol = latent, nrow = ny)
    colnames(d) <- paste0("latent_", seq_len(latent)) 
    new_pfc <- phyf::pf_kronecker(d, pfc)
    # y_latent <- tibble_block(as.data.frame(as.matrix(d)),
    #                          list(dplyr::tibble(y = rep(NA, length(pfc)))))
    new_pfc <- dplyr::tibble("copy_latent_pfc_{num}" := new_pfc)
    # y_latent <- dplyr::bind_cols(y_latent, dplyr::tibble(latent_pfc = new_pfc))
  } else {
    new_pfc <- dplyr::tibble()
  }
  
  new_pfc
  
}

make_latent_y <- function(latent, ylen) {
  
  if(latent > 0) {
    d <- Matrix::.sparseDiagonal(latent, shape = "g")
    colnames(d) <- paste0("latent_", seq_len(latent)) 
    # new_pfc <- phyf::pf_kronecker(d, pfc)
    y_latent <- tibble_block(as.data.frame(as.matrix(d)),
                             list(dplyr::tibble(y = rep(NA, ylen))))
    
  } else {
    y_latent <- dplyr::tibble()
  }
  
  y_latent
  
}

make_y_data <- function(pfc, latent, ny, ylen, num) {
  d <- Matrix::.sparseDiagonal(ny, shape = "g")
  colnames(d) <- paste0("y_", seq_len(ny))
  if(latent == 0) {
    new_pfc <- phyf::pf_kronecker(d, pfc)
    # y_latent <- tibble_block(as.data.frame(as.matrix(d)),
    #                          list(dplyr::tibble(y = rep(NA, length(pfc)))))
    # y_latent <- dplyr::bind_cols(y_latent, dplyr::tibble(latent_pfc = new_pfc))
    new_pfc <- dplyr::tibble("y_pfc_{num}" := new_pfc)
  } else {
    new_pfc <- dplyr::tibble()
  }
  
  new_pfc
  
}

expand_pfc_to_y <- function(pfc, latent, ny) {
  
  if(latent > 0) {
    n <- latent
  } else {
    n <- ny
  }
  
  do.call(cbind, rep(list(phyf::pf_as_sparse(pfc)), n))
  
}

expand_y_to_y <- function(y, latents) {
  
  new_y <- purrr::imap_dfr(y,
                           ~ dplyr::tibble("{.y}" := .x))
  
  latent <- sum(unlist(latents))
  latent_y <- purrr::map_dfr(paste0("latent_", seq_len(latent)),
                              ~ dplyr::tibble("{.x}" := rep(NA, nrow(y))))
  
  dplyr::bind_rows(new_y, latent_y)
  
  
  
}

expand_A_to_y <- function(phy_A, x_A, ny, latents) {
  
  latent <- unlist(latents) > 0
  nlatent <- sum(unlist(latents))
  
  # phy_A_y <- purrr::map_if(phy_A, latent,
  #                          A_zero_out)
  
  phy_A_y <- phy_A
  
  phy_A_latent <- purrr::map_if(phy_A, !latent,
                           A_zero_out)
  
  x_A_y <- x_A
  x_A_latent <- A_zero_out(x_A)
  
  new_A_y <- do.call(cbind, c(list(x_A_y), phy_A_y))
  new_A_latent <- do.call(cbind, c(list(x_A_latent), phy_A_latent))
  
  new_A_y <- do.call(rbind, rep(list(new_A_y), ny))
  new_A_latent <- do.call(rbind, rep(list(new_A_latent), nlatent))
  
  rbind(new_A_y, new_A_latent)
  
  # 
  #   
  # 
  # 
  # A <- cbind(phy_A, x_A)
  # 
  # new_A <- do.call(rbind, rep(list(A), ny))
  # 
  # 
  # 
  # if(latent == 0) {
  #   return(list(x = dat, y = new_y, A = new_A))
  # }
  # 
  # latent_y <- purrr::map_dfr(paste0("latent_", seq_len(latent)),
  #                             ~ dplyr::tibble("{.x}" := rep(NA, nrow(y))))
  # 
  # new_new_y <- dplyr::bind_rows(new_y, latent_y)
  # 
  # latent_x_A <- x_A
  # latent_x_A[] <- 0
  

  
}

check_inla_dims <- function(y, dat, A, ny, nlatent, npredictors, npfc) {
  
  y_match_rows <- dim(y)[1] == dim(A)[1]
  y_match_cols <- dim(y)[2] == ny + nlatent
  
  dat_match_rows <- dim(dat)[1] == dim(A)[2]
  dat_match_cols <- dim(dat)[2] == npredictors + npfc + nlatent
  
  if(any(c(!y_match_rows, !y_match_cols, !dat_match_rows, !dat_match_cols))) {
    message <- paste("Data dimensions don't match:",
                     ifelse(!y_match_rows, "y has incorrect number of rows.", ""),
                     ifelse(!y_match_cols, "y has incorrect number of columns.", ""),
                     ifelse(!dat_match_rows, "x has incorrect number of rows.", ""),
                     ifelse(!dat_match_cols, "x has incorrect number of columns.", ""),
                     collapse = "\n")
    rlang::abort(message)
  }
  
  invisible(TRUE)
  
}

A_zero_out <- function(A) {
  A[] <- 0
  A
}