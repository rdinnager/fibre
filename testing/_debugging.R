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