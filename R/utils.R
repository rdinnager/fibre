assert_inla <- function() {
  if(!requireNamespace("INLA", quietly = TRUE)) {
    stop("Sorry, the fibre package requires the INLA package, which is not installed on your system. Please visit https://www.r-inla.org/download-install and follow the instructions (INLA is not on CRAN currently)")
  }
}

get_vars <- function(form, envir = environment(form)) {
  if(is.list(envir)) {
    envir <- list2env(envir)
  }

  as.data.frame(mget(all.vars(form), envir = envir,
                            ifnotfound = list(NULL)))
}

check_data_dims <- function(y, dat, datas) {

  equal <- TRUE

  y_n <- nrow(y)

  if(!is.null(dat)) {
    dat_n <- nrow(dat)
  } else {
    dat_n <- y_n
  }

  datas_n <- sapply(datas, function(x) nrow(x$rtp_mat))

  if(any(sapply(datas_n, is.null))) {
    others <- sapply(datas[sapply(datas_n, is.null)], function(x) length(x$rtp_mat))
    if(all(others)) {
      datas_n <- datas_n[!others]
      equal <- TRUE
    } else {
      equal <- FALSE
    }
  }

  if(equal) {
    equal <- length(unique(c(y_n, dat_n, datas_n))) == 1
  }

  if(!equal) {
    stop("There is a mismatch of dimension in the model. Implied number of rows are response: ",
         y_n,
         "; fixed predictors: ",
         dat_n,
         "; phylogenetic random effect(s) (in order): ",
         paste(datas_n, collapse = ", "))
  }

  invisible(NULL)
}

#' @importFrom zeallot %<-%
tibble_block <- function(blocks, tibbles, glue_names = TRUE, add_rownames = TRUE) {

  if(!is.list(tibbles)) {
   tibbles <- list(tibbles)
  }

  tibs <- vctrs::vec_recycle_common(!!!tibbles)
  missing <- purrr::map(blocks,
                        ~ tibs[[.x[.x > 0][1]]] %>%
                          dplyr::mutate(dplyr::across(.fns = na_fill)))

  names(missing) <- colnames(blocks)
  new_tib <- blocks %>%
    dplyr::mutate(dplyr::across(.fns = ~ ifelse(.x == 0, missing[dplyr::cur_column()], tibs[.x])))
  if(add_rownames) {
    new_tib <- new_tib %>%
      dplyr::mutate(.rownames = rownames(blocks))
  }

  names_sep <- if(glue_names) ":" else NULL
  unn <- tidyr::unnest(new_tib, dplyr::everything(), names_sep = names_sep)
  unn
}

na_fill <- function(x) {
  x[] <- NA
  x
}

expand_empty <- function(...) {
  obs <- rlang::list2(...)
  empt <- purrr::map_lgl(obs, vctrs::vec_is_empty)
  obs[empt] <- vctrs::vec_init(obs[empt])
  obs
}

empty_sparse <- function(nrow = 0, ncol = 0) {
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  out <- new("dgCMatrix")
  out@Dim <- as.integer(c(nrow, ncol))
  out@p <- integer(ncol + 1L)
  out
}

spark_hist_with_padding <- function(marginals, n_bins = 16) {

  samps <- purrr::map(marginals,
                      ~INLA::inla.rmarginal(.x, n = 400))
  #rngs <- purrr::map(marginals, ~ range(.x[ , "x"]))
  mins <- purrr::map_dbl(samps, ~min(.x))
  maxs <- purrr::map_dbl(samps, ~max(.x))
  minmax_iv <- ivs::iv(mins, maxs)
  breaks <- seq(min(mins), max(maxs), length.out = n_bins)
  bins <- ivs::iv(breaks[-length(breaks)], breaks[-1])
  overlaps <- ivs::iv_locate_overlaps(minmax_iv, bins)
  covers_bins <- overlaps %>%
    dplyr::group_by(needles) %>%
    dplyr::summarise(count = length(haystack),
                     pad_front = min(haystack) - 1,
                     pad_back = n_bins - max(haystack))
  covers_bins <- covers_bins %>%
    dplyr::mutate(count = ifelse(count == 1, 2, count))
  glyphs <- purrr::map2_chr(samps, covers_bins$count,
                            ~ skimr::inline_hist(.x,
                                                 .y))
  padded <- purrr::pmap_chr(list(covers_bins$pad_front, glyphs, covers_bins$pad_back),
                            ~ paste(c(rep(" ", ..1),
                                      ..2,
                                      rep(" ", ..3)), collapse = ""))

  padded

}

spark_dotplot <- function(coef, n_bins = 16) {
  breaks <- seq(min(c(coef, 0)) - 0.01, max(c(coef, 0)) + 0.01, length.out = n_bins)
  bins <- ivs::iv(breaks[-length(breaks)], breaks[-1])
  overlaps <- ivs::iv_locate_between(coef, bins)
  pnt_sym <- cli::format_inline("{cli::symbol$circle_filled}")
  half <- (nchar(pnt_sym) - 1) / 2
  zero_overlap <- ivs::iv_locate_between(0, bins)
  zero_sym <- "|"
  string <- rep(paste(rep(" ", length(bins) + nchar(pnt_sym) - 1), collapse = ""),
                length(coef))
  substr(string, zero_overlap$haystack[1], zero_overlap$haystack[1]) <- zero_sym
  substr(string, overlaps$haystack - half, overlaps$haystack + half) <- pnt_sym
  string
}

get_families <- function(family, y_names, family_hyper = NULL) {


  if(is.null(family_hyper)) {
    family_hyper <- rep(list(list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))))), length(family))
  } else {
    if(length(family) != length(family_hyper)) {
      rlang::abort("You provided hyper-parameters for the family but it's length does not match the length of family.")
    }
  }

  latent <- grep("latent_", y_names)
  non_latent <- setdiff(seq_along(y_names), latent)
  if(length(family) == length(non_latent)) {
    families <- c(family, rep("gaussian", length(latent)))
    family_hyper <- c(family_hyper, rep(list(list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))))), length(latent)))
  } else {
    if(length(family) == 1) {
      families <- c(rep(family, length(non_latent)),
                    rep("gaussian", length(latent)))
      family_hyper <- c(rep(family_hyper, length(non_latent)),
                        rep(list(list(hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1))))), length(latent)))
    } else {
      rlang::abort("family has incorrect length. It should have either length 1 or length equal to the number of outcomes.")
    }
  }

  if(length(latent) > 0) {
    family_hyper[seq_along(non_latent)] <- rep(list(list(hyper = list(prec = list(initial = 10, fixed = TRUE)))), length(non_latent))
  }

  if(family == "binomial") {
    families <- "binomial"
    family_hyper <- list()
  }

  list(family = families, hyper = family_hyper)
}
