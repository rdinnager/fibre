#
Q_criterion <- function(residuals, x, e, n_y) {
  res <- Matrix::t(do.call(rbind, replicate(n_y, cbind(1, x)))) %*% (residuals * e)
  res <- abs(2 * (res / (nrow(x) * n_y)))
  maxes <- apply(res, 2, max)
  maxes
}

resid_sig <- function(res, n_y, phy) {

  Rphylopars::fast.SSC(matrix(res, ncol = n_y), tree = phy,
                       niter = 10)

}

glmnet_metrics <- function(fit, x, y, phy, alpha = 1, penalty.factor = NULL,
                           ncores, verbose, what = c("phylosig",
                                                     "eff_noise"),
                           n_y = 1){

  what <- match.arg(what, c("phylosig",
                            "eff_noise"), several.ok = TRUE)

  if(verbose > 1) {
    cli::cli_progress_message("Calculating metrics along glmnet solution path...")
  }

  if(!is.null(ncores)) {
    rlang::check_installed("furrr")
    oplan <- future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(oplan))
  }

  lambda <- fit$lambda
  #n_y <- ncol(y)
  coef <- coef(fit)
  p <- ncol(x) #* n_y

  if(is.list(coef)) {
    yhat <- lapply(coef, function(b) cbind(1, x) %*% b)
    yhat <- do.call(rbind, yhat)
    y <- as.vector(y)
    #p <- p * length(coef)
    #n_y <- length(coef)
    #print(dim(yhat))
    #print(dim(y))
  } else {
    yhat <- cbind(1, x) %*% coef
    #n_y <- 1
  }

  #df <- fit$df#*n_y
  n <- fit$nobs
  #nlambda <- length(lambda)

  #loocv <- numeric(nlambda)
  #plot(penalty.factor)
  if(!is.null(penalty.factor)) {
    #scal <- penalty.factor * (n / sum(penalty.factor)) #* p
    scal <- (penalty.factor / (sum(penalty.factor))) * p * n
  } else {
    scal <- 1
  }

  # if(alpha == 0){
  #   xs <- x
  #
  #   I <- diag(ncol(x))
  #   xx <- t(xs)%*%xs
  #   for(i in 1:nlambda){
  #
  #     aux <- solve(xx + I * lambda[i] * scal * n)
  #     hatm <- xs%*%aux%*%t(xs)
  #     df[i] <- sum(diag(hatm))
  #
  #     ymat <- matrix(y, ncol = 1)
  #     onemhat <- diag(n) - hatm
  #
  #     loocv[i] <- t(ymat) %*% onemhat %*% (diag(n) * (diag(onemhat)^-2)) %*% onemhat %*% ymat
  #
  #
  #
  #     # B <- diag(1 / (1 - diag(hatm)))
  #     # ff <- (B%*%(diag(ncol(hatm)) - hatm))%*%matrix(y, ncol = 1)
  #     # loocv[i] <- sum(ff^2) / n
  #
  #
  #
  #   # df2 <- df
  #   # t2 <- system.time({
  #   # svd <- sparsesvd(x)
  #   # D <- svd$d^2
  #   # for(i in 1:length(lambda)){
  #   #   aux <- sum(D * (1 / (D + lambda[i] * scal)))
  #   #   df2[i] <- sum(diag(xs%*%aux%*%t(xs)))
  #   # }
  #   # })
  #   }
  #
  #
  # }


  residuals <- (as.vector(y) - yhat)

  if("phylosig" %in% what) {

    resid_sig <- furrr::future_map(purrr::array_branch(as.matrix(residuals), 2),
                                   ~ resid_sig(.x, n_y = n_y, phy = phy),
                                   .progress = verbose > 1,
                                   .options = furrr::furrr_options(seed = TRUE))
    #apply(residuals, 2, resid_sig, n_y = n_y, phy = phy)
    resid_p <- purrr::map_dbl(resid_sig, "pvalue") #sapply(resid_sig, function(x) x$pvalue)
    resid_s <- purrr::map_dbl(resid_sig, "scaled.SSC") #sapply(resid_sig, function(x) x$scaled.SSC)
    resid_sc <- purrr::map_dbl(resid_sig, "SSC")
    resid_interval <- which(resid_p > 0.05 & resid_p < 0.95 & resid_s < 1)
    resid_best <- min(which(resid_p > 0.95 & resid_s < 1))

    if(!is.finite(resid_best)) {
      resid_best <- which.min(resid_s)
    }

    cli::cli_progress_message("finished residual analysis...")
  } else {
    resid_interval <- NULL
    resid_best <- NULL
    resid_p <- NULL
    resid_s <- NULL
  }

  if("eff_noise" %in% what) {

    n <- nrow(residuals)
    scaled_x <- Matrix::t(Matrix::t(x) * (scal))
    Q <- furrr::future_map(1:100,
                           ~ Q_criterion(residuals, scaled_x, rnorm(n), n_y),
                           .progress = verbose > 1,
                           .options = furrr::furrr_options(seed = TRUE))
    Q <- do.call(rbind, Q)
    qQ <- purrr::map_dbl(purrr::array_branch(Q, 2), ~ quantile(.x, probs = 1 - 0.5))

    wmin <- which(qQ <= (lambda * n))
    if(length(wmin) > 0) {
      eff_noise <- min(lambda[wmin])
    } else {
      eff_noise <- lambda[length(lambda)]
    }

    message("finished Q analysis...")
  } else {
    eff_noise <- NULL
    qQ <- NULL
  }

  #mse <- Matrix::colMeans(residuals^2, na.rm = TRUE)
  #sse <- Matrix::colSums(residuals^2, na.rm = TRUE)

  # nvar <- df + 1
  # bic <- 2 * n * log(mse) + nvar * log(n)
  # aic <- 2 * n * log(mse) + 2*nvar
  # aicc <- aic + (2 * nvar * (nvar+1)) / (n - nvar - 1)
  #hqc <- n * log(mse) + 2 * nvar * log(log(n))

  #sigma <- sqrt(sse / max(n-df, 1))

  #  tLL <- fit$nulldev - fit$nulldev * (1 - fit$dev.ratio)
  # k <- df
  # aicc <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
  # aic <- -tLL + 2 * k
  # bic <- log(n * n_y) * k - tLL
  #
  # bic_l <- bic + 2 * lchoose(p, k)

  # k <- df2 + 1
  # aicc2 <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
  # aic2 <- -tLL + 2 * k
  # bic2 <- log(n) * k - tLL

  #rate_est <- (sigma^2) / (lambda * n)

  list(eff_noise = eff_noise, eff_noise_vec = qQ,
       phylosig_interval = resid_interval,
       phylosig_best = resid_best,
       phylosig_vec_p = resid_p, phylosig_vec_ssc = resid_s,
       phylosig_sc = resid_sc,
       lambda_adj = lambda * n)

  # list(bic = bic, aic = aic, aicc = aicc, bic_l = bic_l, sigma = sigma, df = df,
  #      mse = mse, rate_est = rate_est, loocv = loocv, LL = tLL, eff_noise_vec = qQ,
  #      lambda = lambda, eff_noise = eff_noise, resid_interval = resid_interval,
  #      resid_best = resid_best, resid_p = resid_p, resid_sss = resid_s)
}
