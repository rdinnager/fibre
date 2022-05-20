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

prior <- makemyprior_gui(prior)

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


