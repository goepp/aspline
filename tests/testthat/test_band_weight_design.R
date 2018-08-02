context("band_weight_design")

library(tidyverse)
library(splines2)
library(boot)
library(lattice)
load_all()

data(coal)
coal_year <- coal %>%
  mutate(year = floor(date))
set.seed(0)
count <- coal_year %>%
  count(year) %>%
  mutate(year = as.integer(year))
x <- count$year[seq(1, nrow(count), 2)]
y <- count$n[seq(1, nrow(count), 2)]
n <- length(x)
knots <- seq(x[1] + diff(x)[1] / 2,
             x[n] - diff(x)[1] / 2,
             9)
degree <- 1
pen <- 0
X <- bSpline(x, knots = knots, intercept = TRUE, degree = degree)
levelplot(X == 0)
plot(x, y, col = "gray", ylim = c(0, 7))
abline(v = knots, col = "darkgreen")

# pen <- 10 ^ seq(-3, 3, length = 50)
comp <- block_design(X, degree)
B <- comp$B
alpha <- comp$alpha

set.seed(0)
par <- rnorm(ncol(X))
w <- rep(1, ncol(X) - degree - 1)
g_inv <- exp
g_p <- function(x) return(1 / x)

## hessian_solver_glm
Omega <- diag(as.vector(X %*% par))
XWX <- t(X) %*% Omega %*% X
XWX_band <- weight_design_band(as.vector(X %*% par), alpha, B)
nrow(XWX) == nrow(XWX_band)
norm(mat2rot(XWX) - XWX_band)
all.equal(mat2rot(XWX), XWX_band)
## weight_design_band seems erroneous: the dimensions do not agree

## Diagnosis of weight_design
comp <- block_design(X, degree)
B <- comp$B
alpha <- comp$alpha

levelplot(X == 0)
levelplot(t(X) %*% X == 0)

first_non_zero <- function(vect) {
  return(which(vect != 0)[1])
}
shift_level <- apply(X, 1, first_non_zero) %>% unname()
shift_level

add_epsi_old <- function(ind_row, X, degree) {
  X_epsi <- X[ind_row, ]
  shift_level <- apply(X, 1, first_non_zero) %>% unname()
  m <- degree + 1
  index_non_zero <- which(X[ind_row, ] != 0)
  if (ind_row == 1) {
    X_epsi[1:m] <- X_epsi[1:m] + 1e-100
  } else if (diff(shift_level)[ind_row] != 0) {
    if (length(index_non_zero) == m) {
      # do nothing
    } else if (length(index_non_zero) > m) {
      stop("error: band of design matrix is larger than degree")
    } else if (length(index_non_zero) < m) {
      X_epsi[index_non_zero[1]:(index_non_zero[1] + m)] <-
        X_epsi[index_non_zero[1]:(index_non_zero[1] + m)] + 1e-100
    }
  }
  X_epsi
}
add_epsi <- function(X, degree) {
  X_epsi <- X
  shift_level <- apply(X, 1, first_non_zero) %>% unname()
  m <- degree + 1
  for (ind_row in 1:nrow(X)) {
    index_non_zero <- which(X[ind_row, ] != 0)
    if (ind_row == 1) {
      X_epsi[ind_row, 1:m] <- X_epsi[ind_row, 1:m] + 1e-100
    } else if (ind_row == nrow(X)) {
      X_epsi[ind_row, (ncol(X) - m + 1):ncol(X)] <-
        X_epsi[ind_row, (ncol(X) - m + 1):ncol(X)] + 1e-100
    } else if (diff(shift_level)[ind_row - 1] != 0) {
      if (length(index_non_zero) == m) {
        # do nothing
      } else if (length(index_non_zero) > m) {
        stop("error: band of design matrix is larger than degree")
      } else if (length(index_non_zero) < m) {
        X_epsi[ind_row, index_non_zero[1]:(index_non_zero[1] + m - 1)] <-
          X_epsi[ind_row, index_non_zero[1]:(index_non_zero[1] + m - 1)] + 1e-100
      }
    }
  }
  X_epsi
}

# X_epsi <- add_epsi_2(X, degree)
# levelplot(X_epsi == 0)
# m <- degree + 1
# X_epsi <- '[<-'(
#   '[<-'(X, 1, 1:m, 1e-100),
#   nrow(X), (ncol(X) - m + 1):ncol(X), 1e-100)
# levelplot(X_epsi == 0)
# support_by_col <- which(X_epsi != 0, arr.ind = TRUE)
# support_by_row <- support_by_col[order(support_by_col[, 1]), ]
# X0 <- X * 0
# X0[support_by_row] <-  1
# decalage <- apply(X0, 1, function(a) which(a != 0)[1])
# alpha <- c(which(c(1, diff(decalage)) == 1), nrow(X) + 1)
# B <- matrix(X[support_by_row], nrow(X), m, byrow = TRUE)
# list(B = B, alpha = alpha)
