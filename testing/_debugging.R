library(phyf)
library(tidyverse)
library(fibre)

data("vert_bmr")

fibre_mod <- fibre(scale(lnBMR) + scale(lnMass) + scale(lnMass2) ~ bre_brownian(phlo,
                                                               standardise = FALSE),
                   data = vert_bmr)

autoplot(fibre_mod)

fit <- fibre(scale(lnBMR) + scale(lnMass) + scale(lnMass2) ~ bre_brownian(phlo,
                                                                          standardise = FALSE),
             data = vert_bmr,
             engine = "glmnet",
             family = "mgaussian",
             verbose = 2,
             ncores = 3)

fit1 <- fibre(scale(lnBMR) + scale(lnMass) + scale(lnMass2) ~ bre_brownian(phlo,
                                                                          standardise = FALSE),
             data = vert_bmr,
             engine = "glmnet",
             family = "gaussian",
             verbose = 2,
             ncores = 3,
             engine_options = list(alpha = 0.0))

fit2 <- fibre(scale(lnBMR) + scale(lnMass) + scale(lnMass2) ~ bre_brownian(phlo,
                                                                          standardise = FALSE),
             data = vert_bmr,
             engine = "glmnet",
             family = "mgaussian",
             verbose = 2,
             engine_options = list(alpha = 0.0))

stacked <- vert_bmr %>%
  mutate(lnBMR_st = scale(lnBMR),
         lnMass_st = scale(lnMass),
         lnMass2_st = scale(lnMass2)) %>%
  select(lnBMR_st, lnMass_st, lnMass2_st, phlo) %>%
  pivot_longer(cols = -phlo, names_to = "var", values_to = "val") %>%
  mutate(fact = pf_as_pfc(var),
         new_phy = pf_row_kron(fact, phlo))

plot(pf_as_phylo(stacked$new_phy), type = "f")
tree_test <- pf_as_phylo(stacked$new_phy)

fit3 <- fibre(val ~ bre_brownian(new_phy,
                                 standardise = FALSE),
             data = stacked,
             engine = "glmnet",
             family = "gaussian",
             verbose = 2,
             engine_options = list(alpha = 0.01,
                                   what = c("eff_noise", "phylosig")))


tt <- pf_mean_edge_features(vert_bmr)
re <- fibre_mod$random$phlo$`0.5quant` / sqrt(tt)
re1 <- fibre_mod$random$phlo$`0.5quant`
re1[abs(re1) < 0.001] <- 0
re2 <- re1 / sqrt(tt)

plot(re ~ fit2$random$phlo$coef)
plot(re2 ~ fit2$random$phlo$coef)
abline(0, 1)
plot(re2 ~ fit1$random$phlo$coef)
abline(0, 1)

cor(re2, fit2$random$phlo$coef)
cor(re2, fit1$random$phlo$coef)

test<-predict(fibre_mod)
test2 <- purrr::map_dfc(test, ".pred_mean")
tips <- (vert_bmr$is_tip)
yhat <- test2 %>% filter(tips) %>% purrr::flatten_dbl()
y <- na.omit(fibre_mod$extras$data$y$y)
plot(yhat, y)

yhat2 <- unlist(fit2$saved_predictions[tips, -1])
plot(yhat2, y)

test <- predict(fibre_mod)
plot(test[[2]]$.pred_mean ~ vert_bmr$lnMass)

data(bird_beak_codes)
mod <- fibre(latent_code_1 + latent_code_2 +
               latent_code_3 ~ bre_brownian(phlo) + X,
             data = bird_beak_codes)


x <- mod$dat
A <- mod$A
ys_all <- mod$y

set.seed(375)

mod <- fibre(latent_code_1 + latent_code_2 +
               latent_code_3 ~ bre_brownian(phlo),
             data = bird_beak_codes,
             engine_options = list(verbose = TRUE))

y <- bird_beak_codes %>%
  select(starts_with("latent_")) %>%
  as.matrix()

mod <- fibre(y ~ bre_brownian(phlo) + X,
             data = bird_beak_codes %>% mutate(y = y))


data.same.len <- INLA::inla.stack.data(mod[[1]])
data.same.len$scale.log.Mass...1.. <- rep(Inf, length(data.same.len[[2]]))
names(data.same.len)[1] <- "y..fake"
new.fix.formula <- y..fake ~ X.Intercept._1 - 1
contrasts <- NULL
inla.na.action <- function(x, ...) {
            for (k in seq_along(x)) {
                if ((is.numeric(x[[k]]) || inla.is.matrix(x[[k]])) && !is.factor(x[[k]])) {
                    x[[k]][is.na(x[[k]])] <- 0
                }
            }
            return(na.pass(x))
        }


test <- MatrixModels::model.Matrix(
  new.fix.formula,
  data = model.frame(new.fix.formula, data.same.len, na.action = inla.na.action),
  contrasts.arg = contrasts, sparse = TRUE
  )
