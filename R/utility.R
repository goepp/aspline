#' Transform a Spline Design Matrix in block compressed form
#'
#' @param X The design matrix, as given by \code{splines2::bSpline}.
#' @param degree Degree of the spline regression, as used in function \code{splines2::bSpline}.
#' @return A matrix \code{B} with all non-zero entries of \code{X} and a vector of indices \code{alpha} representing the positions of the non-zero blocks of \code{X}.
#' @export
block_design <- function(X, degree) {
  m <- degree + 1
  X_epsi <- '[<-'(
    '[<-'(X, 1, 1:m, 1e-100),
    nrow(X), (ncol(X) - m + 1):ncol(X), 1e-100)
  support_by_col <- which(X_epsi != 0, arr.ind = TRUE)
  support_by_row <- support_by_col[order(support_by_col[, 1]), ]
  X0 <- X * 0
  X0[support_by_row] <-  1
  decalage <- apply(X0, 1, function(a) which(a != 0)[1])
  alpha <- c(which(c(1, diff(decalage)) == 1), nrow(X) + 1)
  B <- matrix(X[support_by_row], nrow(X), m, byrow = TRUE)
  list(B = B, alpha = alpha)
}
