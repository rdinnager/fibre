library(INLA)
library(fibre)
library(makemyprior)

library(phytools)

set.seed(230045)
## transition matrix for the discrete trait simulation
Q <- matrix(c(-1, 1, 1, -1), 2, 2)
## simulated tree
tree <- pbtree(n = 100, scale = 1, d = 0.2)

rates <- setNames(c(1, 20), 1:2)

x <- simulate_traits(tree, rate_model = "discrete", temp_trend_rates = 0,
                     rate_change = Q, rates = rates, internal = FALSE)
## add a bit of noise
x <- x + rnorm(length(x), 0, 0.5)

tree2 <- tree
tree2$edge.length <- sqrt(tree2$edge.length)

A_1 <- make_root2tip(tree2, "tips", order = "first")
A_2 <- make_root2tip(tree2, "tips", order = "second")

## scale_these
A_1 <- A_1 / sqrt(typical_variance(A_1 %*% t(A_1)))
A_2 <- A_2 / sqrt(typical_variance(A_2 %*% t(A_2)))

id1 <- 1:ncol(A_1)
id2 <- 1:ncol(A_2)
intercept <- rep(1, length(x))

data_stack <- inla.stack(data = list(x = x),
                          A = list(1, A_1, A_2),
                          effects = list(root = intercept, id1 = id1, id2 = id2))

dat <- inla.stack.data(data_stack)
A <- inla.stack.A(data_stack)

orig_model <- inla(x ~ 0 + root + f(id1, model = "iid") + f(id2, model = "iid"),
                   data = dat,
                   control.predictor = list(A = A, compute = TRUE))

summary(orig_model)

## Now to make an equivalent model at the observation scale (no A matrices)
C1 <- solve(A_1 %*% t(A_1))
C2 <- solve(A_2 %*% t(A_2))

oid1 <- 1:length(x)
oid2 <- 1:length(x)

odat <- list(x = x, id1 = oid1, id2 = oid2, root = intercept)

obs_model <- inla(x ~ 0 + root + f(id1, model = "generic0", Cmatrix = C1) +
                    f(id2, model = "generic0", Cmatrix = C2),
                  data = odat,
                  control.predictor = list(compute = TRUE))

summary(obs_model)

## These appear to be more or less equivalent model as expected
## I chalk up the difference for the one extremely small effect
## as numerical accuracy due to the A to C conversion.
## Note that the marginal log-likelihood is different because
## "generic0" does not add the factor from the Cmatrix determinant
## see how close mlik is after accounting for this:
orig_model$mlik - 0.5*log(det(C1)) - 0.5*log(det(C2))
obs_model$mlik ## it is pretty close

## Now we specify the above model using makemyprior

prior <- make_prior(x ~ 0 + root +
                      mc(id1, "generic0", constr = FALSE, Cmatrix = C1) +
                      mc(id2, "generic0", constr = FALSE, Cmatrix = C2),
                    data = odat,
                    covariate_prior = list(root = c(0, 1000)))

lower_chol <- prior$lower_cholesky

#prior <- makemyprior_gui(prior)

readr::write_rds(prior, "testing/prior.rds")
prior <- readr::read_rds("testing/prior.rds")
prior$lower_cholesky <- lower_chol

obs_model_priors <- inference_inla(prior,
                                control.predictor = list(compute = TRUE))

summary(obs_model_priors$inla)

## This actually looks really nice, the estimate seem good. The iid effect is closer to the truth,
## and the second order effect is no longer shrunk to oblivion. And now we even have a few more
## effective replicates, which is nice (though still not many). The marginal likelihood has improved
## as well!

## Now how do we modify the makemyprior object to fit the A matrix equivalent model?
## I want this so that I can get the estimates of the random effects for the A matrix columns.
## I figure the priors for these variance components can be used as is in the original model,
## don't see anyway to extract the prior, looks like it uses make_jpr internally..
## let's try playing around with internal functions:

jpr_dat <- makemyprior::make_eval_prior_data(prior)
jpr_dat2 <- makemyprior:::make_inla_data_object(prior, list(control.predictor = list(compute = TRUE)))

## those appear to be more or less the same.. let's try this

orig_model2 <- inla(x ~ 0 + root + f(id1, model = "iid") + f(id2, model = "iid"),
                   data = dat,
                   control.predictor = list(A = A, compute = TRUE),
                   control.expert = list(jp = makemyprior:::make_jpr(jpr_dat)))

summary(orig_model2)

## It worked!!! I guess now all we need is for make_jpr to be exported?
## Maybe there is a way to use INLA::inla.jp.define() and makemyprior::eval_joint_prior()?



############## setup a function to do all that on any data #############

library(ape)
library(readr)
library(dplyr)
library(glmnet)
library(Matrix)
library(sparseMatEst)
library(qs)

bt_res <- read_rds("extdata/rateShiftInformation.RDS")
bird_tree <- bt_res$Phylogeny
bird_dat <- read_csv("extdata/AVONET3_BirdTree.csv")

bird_dat <- bird_dat %>%
  mutate(species = gsub(" ", "_", Species3)) %>%
  filter(species %in% bird_tree$tip.label)

tree <- bird_tree
x <- log(bird_dat$Mass)
names(x) <- bird_dat$species


tree2 <- tree
tree2$edge.length <- sqrt(tree2$edge.length)

A_1 <- make_root2tip(tree2, "tips", order = "first")
A_2 <- make_root2tip(tree2, "tips", order = "second")

#A_1 <- A_1[-which(rownames(A_1) == as.character(Ntip(tree) + 1)), ]
#A_2 <- A_2[-which(rownames(A_2) == as.character(Ntip(tree) + 1)), ]

A_1 <- A_1[match(names(x), rownames(A_1)), ]
A_2 <- A_2[match(names(x), rownames(A_2)), ]

## scale_these
A_1 <- A_1 / sqrt(typical_variance(A_1 %*% t(A_1)))
A_2 <- A_2 / sqrt(typical_variance(A_2 %*% t(A_2)))

id1 <- 1:ncol(A_1)
id2 <- 1:ncol(A_2)
intercept <- rep(1, length(x))

C1 <- solve(A_1 %*% t(A_1))
C2 <- solve(A_2 %*% t(A_2))

readr::write_rds(C1, "C1.rds")
readr::write_rds(C1, "C2.rds")

x <- sim$dat
tree <- sim$tree

run_full_mod <- function(x, tree = NULL, C1 = NULL, C2 = NULL) {

  if(is.matrix(x)) {
    nam <- rownames(x)
    x <- as.data.frame(x)
  } else {
    nam <- names(x)
  }

  if(is.null(C1)) {
    tree2 <- tree
    tree2$edge.length <- sqrt(tree2$edge.length)

    A_1 <- make_root2tip(tree2, "tips", order = "first")
    A_2 <- make_root2tip(tree2, "tips", order = "second")

    #A_1 <- A_1[-which(rownames(A_1) == as.character(Ntip(tree) + 1)), ]
    #A_2 <- A_2[-which(rownames(A_2) == as.character(Ntip(tree) + 1)), ]

    A_1 <- A_1[match(nam, rownames(A_1)), ]
    A_2 <- A_2[match(nam, rownames(A_2)), ]

    ## scale_these
    A_1 <- A_1 / sqrt(typical_variance(A_1 %*% t(A_1)))
    A_2 <- A_2 / sqrt(typical_variance(A_2 %*% t(A_2)))

    C1 <- solve(A_1 %*% t(A_1))
    C2 <- solve(A_2 %*% t(A_2))

  } else {
    stop("Either tree or C1 and C2 must be specified")
  }

  id1 <- 1:ncol(C1)
  id2 <- 1:ncol(C2)
  intercept <- rep(1, length(nam))
  oid1 <- 1:length(nam)
  oid2 <- 1:length(nam)


  if(is.data.frame(x)) {
    odat <- list(V1 = x$V1, V2 = x$V2, id1 = oid1, id2 = oid2, root = intercept,
                 C1 = C1, C2 = C2)
    # obs_model <- inla(V2 ~ 0 + root + V1 + f(id1, model = "generic0", Cmatrix = C1) +
    #                     f(id2, model = "generic0", Cmatrix = C2),
    #                   data = odat,
    #                   control.predictor = list(compute = TRUE))
    data_stack <- inla.stack(data = list(V2 = x$V2),
                             A = list(1, A_1, A_2, 1),
                             effects = list(root = intercept, id1 = 1:ncol(A_1), id2 = 1:ncol(A_2), V1 = x$V1))

    dat <- inla.stack.data(data_stack)
    A <- inla.stack.A(data_stack)
  } else {
    odat <- list(x = x, id1 = oid1, id2 = oid2, root = intercept,
                 C1 = C1, C2 = C2)
    # obs_model <- inla(x ~ 0 + root + f(id1, model = "generic0", Cmatrix = C1) +
    #                     f(id2, model = "generic0", Cmatrix = C2),
    #                   data = odat,
    #                   control.predictor = list(compute = TRUE))
    data_stack <- inla.stack(data = list(x = x),
                             A = list(1, A_1, A_2),
                             effects = list(root = intercept, id1 = 1:ncol(A_1), id2 = 1:ncol(A_2)))

    dat <- inla.stack.data(data_stack)
    A <- inla.stack.A(data_stack)

  }

  ## Now we specify the above model using makemyprior

  if(is.data.frame(x)) {
    prior <- make_prior(V2 ~ 0 + root + V1 +
                          mc(id1, "generic0", constr = FALSE, Cmatrix = C1) +
                          mc(id2, "generic0", constr = FALSE, Cmatrix = C2),
                        data = odat,
                        prior = list(tree = "id1_id2_eps = (eps,id1_id2); id1_id2 = (id1,id2)",
                                     V = list(id1_id2_eps = list(prior = "jeffreys")),
                                     w = list(id1_id2_eps = list(prior = "pc0", param = 0.5),
                                              id1_id2 = list(prior = "pc1", param = 0.5))),
                        covariate_prior = list(root = c(0, 1000),
                                               V1 = c(0, 100)))
  } else {
    prior <- make_prior(x ~ 0 + root +
                          mc(id1, "generic0", constr = FALSE, Cmatrix = C1) +
                          mc(id2, "generic0", constr = FALSE, Cmatrix = C2),
                        data = odat,
                        prior = list(tree = "id1_id2_eps = (eps,id1_id2); id1_id2 = (id1,id2)",
                                     V = list(id1_id2_eps = list(prior = "jeffreys")),
                                     w = list(id1_id2_eps = list(prior = "pc0", param = 0.5),
                                              id1_id2 = list(prior = "pc1", param = 0.5))),
                        covariate_prior = list(root = c(0, 1000)))
  }

  # lower_chol <- prior$lower_cholesky
  #
  # prior <- makemyprior_gui(prior)

  # readr::write_rds(prior, "testing/prior.rds")
  # prior <- readr::read_rds("testing/prior.rds")
  # prior$lower_cholesky <- lower_chol

  # obs_model_priors <- inference_inla(prior,
  #                                    control.predictor = list(compute = TRUE))
  #
  # summary(obs_model_priors$inla)

  ## This actually looks really nice, the estimate seem good. The iid effect is closer to the truth,
  ## and the second order effect is no longer shrunk to oblivion. And now we even have a few more
  ## effective replicates, which is nice (though still not many). The marginal likelihood has improved
  ## as well!

  ## Now how do we modify the makemyprior object to fit the A matrix equivalent model?
  ## I want this so that I can get the estimates of the random effects for the A matrix columns.
  ## I figure the priors for these variance components can be used as is in the original model,
  ## don't see anyway to extract the prior, looks like it uses make_jpr internally..
  ## let's try playing around with internal functions:

  jpr_dat <- makemyprior::make_eval_prior_data(prior)
  jpr_dat2 <- makemyprior:::make_inla_data_object(prior, list(control.predictor = list(compute = TRUE)))

  ## those appear to be more or less the same.. let's try this

   if(is.data.frame(x)) {
     orig_model2 <- inla(V2 ~ 0 + root + V1 + f(id1, model = "iid") + f(id2, model = "iid"),
                         data = dat,
                         control.predictor = list(A = A, compute = TRUE),
                         control.expert = list(jp = makemyprior:::make_jpr(jpr_dat)),
                         control.compute = list(waic = TRUE))
   } else {
     orig_model2 <- inla(V2 ~ 0 + root + f(id1, model = "iid") + f(id2, model = "iid"),
                         data = dat,
                         control.predictor = list(A = A, compute = TRUE),
                         control.expert = list(jp = makemyprior:::make_jpr(jpr_dat)),
                         control.compute = list(waic = TRUE))
   }

  #summary(orig_model2)

  orig_model2
}


#### felsenstein's worst case ############

library(bayou)
library(treeplyr)
library(geiger)
library(ape)
library(phytools)
library(selectiveInference)

felsTree <- function(n, stemL = 1, lambda = 1, S = 1.5){
  phy1 <- ape::reorder.phylo(sim.bdtree(b=1, d=0, stop="taxa", n=n), "postorder")
  phy2 <- ape::reorder.phylo(sim.bdtree(b=1, d=0, stop="taxa", n=n), "postorder")
  phy1 <- rescale(phy1, model="lambda", lambda)
  phy2 <- rescale(phy1, model="lambda", lambda)
  phy1$edge.length <- phy1$edge.length/max(branching.times(phy1))
  phy2$edge.length <- phy2$edge.length/max(branching.times(phy2))
  o <- list(phy1=phy1, phy2=phy2)
  phy2$tip.label <- paste("s", (n+1):(2*n), sep="")
  phy <- list(edge=NULL, Nnode=(n-1)*2 + 1, tip.label=c(phy1$tip.label, phy2$tip.label), edge.length=NULL)
  phy1$edge[phy1$edge > n] <- phy1$edge[phy1$edge > n]+(n+1)
  phy2$edge[phy2$edge > n] <- phy2$edge[phy2$edge > n]+2*n
  phy2$edge[phy2$edge <= n] <- phy2$edge[phy2$edge <=n]+n

  phy$edge <- rbind(phy1$edge, phy2$edge, c((2*n+1), (2*n+n+1)), c((2*n+1), (2*n+2)))
  phy$edge.length <- c(phy1$edge.length, phy2$edge.length, rep(stemL, 2))
  class(phy) <- class(o$phy1)
  phy <- ape::reorder.phylo(phy, "postorder")

  phy$edge.length <- phy$edge.length/(max(branching.times(phy)))

  S=10^S
  dat <- sim.char(phy, par=matrix(c(1, 0, 0, 1), ncol=2), model="BM", root=0)[,,1]
  dat2 <- dat

  dev <- MASS::mvrnorm(1, mu=c(0,0), Sigma=matrix(c(S, 0*S, 0*S, S), ncol=2))
  dat2[1:n, ] <- dat2[1:n,] + cbind(rep(dev[1], n), rep(dev[2], n))

  return(list(tree = phy, dat = dat2))
}

test <- felsTree(100)

n_sims <- 2000
size_range <- c(20, 40)

for(i in seq_len(n_sims)) {
  n <- sample.int(size_range[2] - size_range[1], 1) + size_range[1]
  lambda <- runif(1, 0.5, 1)
  s <- rnorm(1, 1.5, 1)

  sim <- felsTree(n, lambda = lambda, S = s)
  res1 <- try(run_full_mod(sim$dat, sim$tree))

  tree2 <- sim$tree
  tree2$edge.length <- sqrt(tree2$edge.length)
  res2 <- try(fibre(V2 ~ V1 + p(tree2), data = sim$dat))
  res3 <- try(inla(V2 ~ V1, data = as.data.frame(sim$dat)))

  pic1 <- pic(sim$dat[,1], sim$tree)
  pic2 <- pic(sim$dat[,2], sim$tree)
  res4 <- try(lm(pic2 ~ pic1 - 1))


  rtp <- make_root2tip(sim$tree)

  edges_nums <- c((Ntip(sim$tree) + 2):(Ntip(sim$tree) + Nnode(sim$tree)), 1:Ntip(sim$tree))
  edge_ord <- match(edges_nums, sim$tree$edge[ , 2])
  lens <- sim$tree$edge.length[edge_ord]

  rtp <- cbind(sim$dat[ , 1], rtp)

  rtp <- scale(rtp, TRUE, c(0.01, lens))

  res5 <- try(glmnet(rtp, sim$dat[ , 2], standardize = FALSE, nlambda = 500, lambda = seq(1, 0, length.out = 500)))
  aics <- try(crits(res5, rtp, sim$dat[ , 2]))
  best_lam <- res5$lambda[which.min(aics$aic2)]
  sig <- aics$sigma[which.min(aics$aic2)]
  n <- nrow(sim$dat)
  beta_hat <- coef(res5, x = rtp, y = sim$dat[ , 2], s = best_lam, exact = TRUE)[-1]
  out <- fixedLassoInf(rtp, sim$dat[ , 2], beta_hat, best_lam * n, sigma = sig, bits = 200)

  qsave(list(res1 = res1, res2 = res2, res3 = res3, res4 = res4, tree = sim$tree, dat = sim$dat,
             n = n, lambda = lambda, s = s),
        paste0(file.path("testing/results", "fels_"), i, ".qs"))

  print(paste("Done", i, "out of", n_sims))
}

crits <- function(fit, x, y, alpha = 1){

  lambda <- fit$lambda
  df <- fit$df
  n <- length(y)

  if(alpha == 0){
    xs = x
    I = diag(ncol(x))
    xx = t(xs)%*%xs
    for(i in 1:length(lambda)){
      aux = solve(xx + I * lambda[i] * n)
      df[i] = sum(diag(xs%*%aux%*%t(xs)))
    }


  }

  coef <- coef(fit)

  yhat <- cbind(1, x) %*% coef
  residuals <- (y - yhat)
  mse <- colMeans(residuals^2)
  sse <- colSums(residuals^2)

  nvar <- df + 1
  bic <- 2 * n * log(mse) + nvar * log(n)
  aic <- 2 * n * log(mse) + 2*nvar
  aicc <- aic + (2 * nvar * (nvar+1)) / (n - nvar - 1)
  #hqc <- n * log(mse) + 2 * nvar * log(log(n))

  sigma <- sqrt(sse / (n-df-1))

  tLL <- fit$nulldev - fit$nulldev * (1 - fit$dev.ratio)
  k <- df
  aicc2 <- -tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
  aic2 <- -tLL + 2 * k
  bic2 <- log(n) * k - tLL

  rate_est <- (sigma^2) / (lambda * n)

  list(bic = bic, aic = aic, aicc = aicc, sigma = sigma, df = df,
       mse = mse, aicc2 = aicc2, bic2 = bic2, aic2 = aic2, rate_est = rate_est)
}


######### analyse results ##############
library(purrr)
library(qs)
library(dplyr)

qs_files <- list.files("testing/results", full.names = TRUE, pattern = ".qs")

grab_data <- function(qs_file) {

  res <- qread(qs_file)
  ci_uncont <- c(res$res3$summary.fixed$`0.025quant`[2], res$res3$summary.fixed$`0.975quant`[2])
  pos_uncont <- all(ci_uncont < 0) | all(ci_uncont > 0)

  ci_simpcont <- c(res$res2$summary.fixed$`0.025quant`[2], res$res2$summary.fixed$`0.975quant`[2])
  pos_simpcont <- all(ci_simpcont < 0) | all(ci_simpcont > 0)

  ci_flexcont <- c(res$res1$summary.fixed$`0.025quant`[2], res$res1$summary.fixed$`0.975quant`[2])
  pos_flexcont <- all(ci_flexcont < 0) | all(ci_flexcont > 0)

  ci_piccont <- confint(res$res4)
  pos_piccont <- all(ci_piccont < 0) | all(ci_piccont > 0)

  rbind(data.frame(mod = "uncont", pos = pos_uncont, lower = ci_uncont[1], upper = ci_uncont[2],
                   n = res$n, l = res$lambda, S = res$s),
        data.frame(mod = "simple", pos = pos_simpcont, lower = ci_simpcont[1], upper = ci_simpcont[2],
                   n = res$n, l = res$lambda, S = res$s),
        data.frame(mod = "flexible", pos = pos_flexcont, lower = ci_flexcont[1], upper = ci_flexcont[2],
                   n = res$n, l = res$lambda, S = res$s),
        data.frame(mod = "pic", pos = pos_piccont, lower = ci_piccont[1], upper = ci_piccont[2],
                   n = res$n, l = res$lambda, S = res$s))

}

dat <- map_dfr(qs_files,
               grab_data)

summ <- dat %>%
  mutate(s_bin = cut(S, breaks = 5)) %>%
  group_by(mod, s_bin) %>%
  summarise(prop_pos = sum(pos) / n())
