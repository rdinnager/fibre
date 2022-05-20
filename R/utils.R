assert_inla <- function() {
  if(!requireNamespace("INLA", quietly = TRUE)) {
    stop("Sorry, the fibre package requires the INLA package, which is not installed on your system. Please visit https://www.r-inla.org/download-install and follow the instructions (INLA is not on CRAN currently)")
  }
}

make_Cmatrix <- function(phy, standardise_var = TRUE, tips_only = FALSE) {

  if(tips_only) {

    A <- MCMCglmm::inverseA(phy, scale = FALSE, nodes = "TIPS")
    Amat <- A$Ainv

  } else {

    tips <- 1:ape::Ntip(phy)
    nodes <- (ape::Ntip(phy) + 1):(ape::Ntip(phy) + ape::Nnode(phy))
    A <- MCMCglmm::inverseA(phy, scale = FALSE)
    Amat <- A$Ainv
    new_names <- rownames(Amat)
    Anodes <- 1:(ape::Nnode(phy) - 1)
    Atips <- seq_along(new_names)[-Anodes]
    new_names[Anodes] <- nodes[as.numeric(gsub("Node", "", new_names[Anodes]))]
    rownames(Amat) <- colnames(Amat) <- new_names

  }

  if(standardise_var) {

    gen_var <- 1 / exp(Matrix::determinant(Amat)$modulus[1]/ape::Ntip(phy))
    Amat <- Amat * gen_var

  }

  Amat
}

make_extraconstr <- function(phy, parents, order) {
  if(!ape::is.binary(phy)) {
      stop("Using a constrained model requires a binary tree.")
    }
    nnode <- ape::Nnode(phy)
    children <- unlist(split(seq_len(nrow(parents)), parents[ , 1]))
    i <- rep(1:nnode, each = 2)
    vals <- rep(c(1, 1), nnode)

    constr_A <- Matrix::sparseMatrix(i = i,
                                    j = children,
                                    dims = c(nnode, nnode + ape::Ntip(phy)),
                                    x = vals)

    constr_A <- constr_A[ , -1]

    if(order == "first") {
      diag(constr_A) <- -2
    }

    #test <- rtp_mat[[1]] %*% t(constr_A)
    constr_e <- rep(0, nrow(constr_A))

    list(A = as.matrix(constr_A), e = constr_e)
}

parse_formula <- function(form, data = NULL) {

  if(is.null(data)) {
    data <- mget(all.vars(form), envir = environment(form),
                 inherits = TRUE)
  } else {
    if(is.matrix(data)) {
      if(!is.null(rownames(data))) {
        node_names <- rownames(data)
      } else {
        node_names <- NA
      }
      data <- as.data.frame(data)
      data$none_names <- node_names
    }

    if(is.list(data)) {
      data <- process_data(data)
    }

    if(any(!all.vars(form) %in% names(data))) {
      extra <- mget(all.vars(form)[!all.vars(form) %in% names(data)],
                    envir = environment(form),
                    ifnotfound = list(NULL),
                    inherits = TRUE)
      extra <- extra[sapply(extra, function(x) inherits(x, "phylo"))]
      data <- c(data, extra)
    }
  }

  y <- get_vars(update.formula(form, . ~ 0), data)
  num_y <- ncol(y)

  ts <- terms.formula(form, specials = c("p", "f", "root", "age"), keep.order = TRUE)
  tls <- rownames(attr(ts, "factors"))
  phybres <- tls[attr(ts, "specials")$p]
  #fs <- tls[attr(ts, "specials")$f]
  if(!is.null(attr(ts, "special")$root)) {
    root_remove <- which(attr(ts, "term.labels") %in% tls[is.null(attr(ts, "special")$root)])
    root <- TRUE
  } else {
    root_remove <- NULL
    root <- FALSE
  }

  if(length(phybres) > 0) {
    .fibre_env$multiple_order <- TRUE
    datas <- lapply(phybres, function(x) eval(parse(text = x), envir = c(data, p = p)))
    .fibre_env$multiple_order <- FALSE
  }
  to_remove <- c(which(attr(ts, "term.labels") %in% phybres),
                 #which(attr(ts, "term.labels") %in% fs),
                 root_remove,
                 NULL)

  if(length(to_remove) == length(attr(ts, "term.labels"))) {
    if(root || attr(ts, "intercept") == 1) {
      new_form <- ~ root
    } else {
      new_form <- ~ 0
    }
  } else {
    new_form <- drop.terms(ts, to_remove, keep.response = TRUE)
    new_form <- formula(delete.response(new_form))
    if(root || attr(ts, "intercept") == 1) {
      new_form <- update(new_form, ~ root + .)
    } else {
      new_form <- update(new_form, ~ 0 + .)
    }
  }

  if(length(setdiff(all.vars(new_form), "root")) > 0) {
    if(root || attr(ts, "intercept") == 1) {
      data$root <- 1
    }
    dat <- get_vars(new_form, data)
  } else {
    if(root || attr(ts, "intercept") == 1) {
      dat <- data.frame(root = rep(1, nrow(y)))
    } else {
      dat <- NULL
    }
  }

  check_data_dims(y, dat, datas)

  data_stack <- make_inla_stack(y, datas, dat)

  parms <- data_stack$names
  inla_forms <- mapply(generate_inla_formula,
                       parms$effect_names,
                       parms$rate_models,
                       parms$Cmats,
                       parms$hypers,
                       parms$constrs)

  resp <- all.vars(form[[2]])

  final_form <- reformulate(c(all.vars(new_form[[2]]),
                              sapply(inla_forms, function(x) Reduce(paste, deparse(x)))),
                            resp,
                            intercept = FALSE,
                            env = emptyenv())

  final_form <- Reduce(paste, deparse(final_form))

  res <- list(formula = final_form, data = data_stack)
  attr(res, "debug") <- datas

  res

}


process_data <- function(data) {

  mats <- sapply(data, is.matrix)
  data[mats] <- lapply(data[mats], as.data.frame)
  dfs <- sapply(data, is.data.frame)
  new_data <- data[!dfs]
  dfs <- unlist(unname(data[dfs]), recursive = FALSE)

  c(new_data, dfs)

}

make_inla_stack <- function(y, datas, dat) {

  effect_id <- sapply(datas, function(x) x$name)
  if(any(duplicated(effect_id))) {
    effect_id <- make.unique(effect_id)
    for(i in seq_along(effects)) {
      names(effects[[i]]$data) <- paste0("phy_", effect_id)
      effects[[i]]$name <- effect_id
    }
  }

  effects <- lapply(datas, function(x) x$data)
  effects <- c(effects, list(x = dat))

  As <- lapply(datas, function(x) x$rtp_mat)
  As <- c(As, list(1))

  tces <- lapply(datas, function(x) x$tces_stack)
  aces <- lapply(datas, function(x) x$aces_stack)

  for(i in seq_along(effect_id)) {
    if(!is.null(tces[[i]])) {
      new_data <- tces[[i]]$data$data
      colnames(new_data) <- colnames(y)
      tces[[i]] <- INLA::inla.stack(data = new_data,
                                    A = tces[[i]]$A,
                                    effects = tces[[i]]$effects$data,
                                    tag = "tces")
      # colnames(tces[[i]]$data$data) <- colnames(y)
    }
    if(!is.null(aces[[i]])) {
      new_data <- aces[[i]]$data$data
      colnames(new_data) <- colnames(y)
      aces[[i]] <- INLA::inla.stack(data = new_data,
                                    A = aces[[i]]$A,
                                    effects = aces[[i]]$effects$data,
                                    tag = "aces")
    }
  }

  Cmats <- lapply(datas, function(x) x$Cmat)
  names(Cmats) <- paste0("Cmat_", effect_id)

  constrs <- lapply(datas, function(x) x$extra_constr)
  names(constrs) <- paste0("constrs_", effect_id)

  hypers <- lapply(datas, function(x) x$hyper)
  names(hypers) <- paste0("hyper_", effect_id)

  rate_models <- sapply(datas, function(x) x$rate_model)

  stack <- INLA::inla.stack(data = list(y = y),
                            A = As,
                            effects = effects,
                            tag = "main")

  tces_stack <- do.call(INLA::inla.stack, tces)
  aces_stack <- do.call(INLA::inla.stack, aces)

  stack <- INLA::inla.stack(stack, tces_stack, aces_stack)
  stack_data <- do.call(INLA::inla.stack.data, c(list(stack), Cmats, constrs, hypers))
  stack_A <- INLA::inla.stack.A(stack)

  list(stack = stack, data = stack_data, A = stack_A, names = list(effect_names = paste0("phy_", effect_id),
                                                                   rate_models = rate_models,
                                                                   Cmats = names(Cmats),
                                                                   constrs = names(constrs),
                                                                   hypers = names(hypers)))

}

generate_inla_formula <- function(effect_name, rate_model,
                                  Cmat, hyper, extraconstr) {

  if(rate_model == "fixed") {
    bquote(.(as.name(effect_name)))
  } else {
    bquote(f(.(as.name(effect_name)),
             model = .(rate_model),
             Cmatrix = .(as.name(Cmat)),
             constr = FALSE,
             hyper = .(as.name(hyper)),
             extraconstr = .(as.name(extraconstr))))
  }
}

make_name <- function(rate_model, rate_order, weights, rate_weights) {
  hasw <- if(!is.null(weights)) "_wghts" else ""
  hasrw <- if(!is.null(rate_weights)) "_rtwghts" else ""
  mod_abbr <- switch(rate_model,
                     iid = "iid",
                     phylogenetic = "phylo",
                     fixed = "fxd")
  ord_abbr <- switch(rate_order,
                     first = "_frst",
                     second = "_scnd")
  name <- paste0("phy_", mod_abbr, ord_abbr, hasw, hasrw)
  name
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
