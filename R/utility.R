#' Transform a Spline Design Matrix in block compressed form
#'
#' @param X The design matrix, as given by \code{splines2::bSpline}.
#' @param degree Degree of the spline regression, as used in function \code{splines2::bSpline}.
#' @return A matrix \code{B} with all non-zero entries of \code{X} and a vector of indices \code{alpha} representing the positions of the non-zero blocks of \code{X}.
block_design <- function(X, degree) {
  m <- degree + 1
  X_epsi <- add_epsi(X, degree)
  support_by_col <- which(X_epsi != 0, arr.ind = TRUE)
  support_by_row <- support_by_col[order(support_by_col[, 1]), ]
  X0 <- X * 0
  X0[support_by_row] <-  1
  decalage <- apply(X0, 1, function(a) which(a != 0)[1])
  alpha <- c(which(c(1, diff(decalage)) == 1), nrow(X) + 1)
  B <- matrix(X[support_by_row], nrow(X), m, byrow = TRUE)
  list(B = B, alpha = alpha)
}
first_non_zero <- function(vect) {
  return(which(vect != 0)[1])
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
