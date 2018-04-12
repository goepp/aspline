#' @useDynLib aspline
#' @importFrom Rcpp sourceCpp
NULL
#'
#' Inverse the hessian and multiply it by the score
#'
#' @param par The parameter vector
#' @param XX_band The matrix \code{t(X) %*% X} where \code{X} is the design matrix. This argument is given
#' in the form of a band matrix, i.e., successive columns represent superdiagonals.
#' @param Xy The vector of currently estimated points \code{t(X) %*% y}, where \code{y} is the y-coordinate of the data.
#' @param pen Positive penalty constant.
#' @param w Vector of weights. Has to be of length
#' @param diff The order of the differences of the parameter. Equals \code{degree + 1} in adaptive spline regression.
#' @return The solution of the linear system
hessian_solver <- function(par, XX_band, Xy, pen, w, diff) {
  if (ncol(XX_band) != diff + 1) stop("Error: XX_band must have diff + 1 columns")
  bandsolve(XX_band + pen * band_weight(w, diff), Xy) - par
}
wridge_solver_fast_2 <- function(XX_band, Xy, degree, pen,
                                 w = rep(1, nrow(XX_band) - degree - 1),
                                 old_par = rep(1, nrow(XX_band)),
                                 maxiter = 1000) {
  for (iter in 1:maxiter) {
    par <- old_par + hessian_solver_fast_2(old_par, XX_band, Xy,
                                           pen, w, diff = degree + 1)
    idx <- old_par != 0
    rel_error <- max(abs(par - old_par)[idx] / abs(old_par)[idx])
    if (rel_error < 1e-8) break
    old_par <- par
  }
  par
}
