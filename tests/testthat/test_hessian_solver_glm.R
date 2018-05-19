context("hessian_solver_glm")

library(tidyverse)
library(splines2)
library(boot)
library(lattice)
load_all()

data(coal)

x <- count$year
y <- count$n
n <- length(x)
# knots <- x[-1] - diff(x)[1] / 2
knots <- seq(x[1] + diff(x)[1] / 2,
             x[n] - diff(x)[1] / 2,
             4)
degree <- 1
pen <- 0
X <- bSpline(x, knots = knots, intercept = TRUE, degree = degree)

plot(x, y, col = "gray", ylim = c(0, 7))
abline(v = knots, col = "darkgreen")

# pen <- 10 ^ seq(-3, 3, length = 50)
comp <- block_design(X, degree)
B <- comp$B
alpha <- comp$alpha

set.seed(0)
par <- rnorm(ncol(X))
w = rep(1, ncol(X) - degree - 1)
g_inv <- exp
g_p <- function(x) 1 / x

## hessian_solver_glm
W <- diag(as.vector(X %*% par))
D <- diff(diag(ncol(X)), differences = degree + 1)
mat_1 <- t(X) %*% W %*% X  + pen * t(D) %*% diag(w) %*% D

## hessian_solver_glm_band
glm_weight <- as.vector(X %*% par)
XWX_band <- cbind(weight_design_band(glm_weight, alpha, B), 0)
mat_2 <- XWX_band + pen * band_weight(w, degree + 1)

# They are the same
cbind(bandsolve::mat2rot(mat_1), 0) - mat_2

## hessian_solver_glm
vect_1 <- t(X) %*% W %*% (y - X %*% par) + t(X) %*% W %*% X %*% par
eta <- X %*% old_par
mu <- exp(eta)
vect <- t(X) %*% W %*% (eta + (y - mu) * g_p(mu))

## hessian_solver_glm_band
vect_2 <- sweep(t(X), MARGIN = 2, glm_weight, `*`) %*% y

# They are the same
vect_1 - vect_2
