library(splines2)
library(microbenchmark)
library(aspline)
library(tidyverse)

set.seed(0)
n <- 100
lower_x <- -1
upper_x <- 3
sigma <- 0.1
x <- runif(n, lower_x, upper_x) %>% sort()
eta <- (sin(x) * cos(x) ^ 2 / x) + (1 - x ^ 2 / 1000) + rnorm(n, sd = sigma)
mu <- exp(eta)
y <- rpois(length(mu), mu)
barplot(y)

knots <- seq(-0.5, 2.5, length = 10)
# knots <- seq(-0.5, 2.5, length = 3)
degree <- 1
pen <- 0
X <- bSpline(x, knots = knots, intercept = TRUE, degree = degree)
comp <- block_design(X, degree)
B <- comp$B
alpha <- comp$alpha
w <- rnorm(ncol(X) - degree - 1)
par <- rnorm(ncol(X))

test_that("Band and Regular functions yield the same result", {
  expect_equal(hessian_solver_glm_band(par, X, y, B, alpha, pen, w, degree, family = family),
               hessian_solver_glm(par, X, y, degree, pen, w, family = family))
})

microbenchmark(hessian_solver_glm_band(par, X, y, B, alpha, pen, w, degree, family = family),
             hessian_solver_glm(par, X, y, degree, pen, w, family = family))

