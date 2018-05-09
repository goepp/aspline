context("band_design.R")

library(splines2)
library(lattice)
library(bandsolve)

set.seed(0)
n <- 14
lower_x <- -1
upper_x <- 3
sigma <- 0.1
x <- runif(n, lower_x, upper_x) %>% sort()
y <- (sin(x) * cos(x) ^ 2 / x) + (1 - x ^ 2 / 1000) + rnorm(n, sd = sigma)
K <- 6
knots <- seq(-1, 3, length = K + 2)[-c(1, K + 2)]
degree <- 1
X <- bSpline(x, knots = knots, intercept = TRUE, degree = degree)

w <- rep(1, nrow(X))
w <- rnorm(nrow(X))

comp <- block_design(X, degree)
bandsolve::mat2rot(t(X) %*% sweep(X, MARGIN = 1, w, `*`)) - weight_design_band(w, comp$alpha, comp$B)


solve(t(X) %*% sweep(X, MARGIN = 1, w, `*`)) - bandsolve(weight_design_band(w, comp$alpha, comp$B))
