#' @useDynLib aspline
#' @importFrom Rcpp sourceCpp
NULL
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
#'
#' Inverse the hessian and multiply it by the score
#'
#' @param par The parameter vector
#' @param XX_band The matrix \eqn{X^T X} where \code{X} is the design matrix. This argument is given
#' in the form of a band matrix, i.e., successive columns represent superdiagonals.
#' @param Xy The vector of currently estimated points \eqn{X^T y}, where \eqn{y} is the y-coordinate of the data.
#' @param pen Positive penalty constant.
#' @param w Vector of weights. Has to be of length
#' @param diff The order of the differences of the parameter. Equals \code{degree + 1} in adaptive spline regression.
#' @return The solution of the linear system: \deqn{(X^T X + pen D^T diag(w) D) ^ {-1} X^T y - par}
#' @export
hessian_solver <- function(par, XX_band, Xy, pen, w, diff) {
  if (ncol(XX_band) != diff + 1) stop("Error: XX_band must have diff + 1 columns")
  bandsolve::bandsolve(XX_band + pen * band_weight(w, diff), Xy) - par
}
#' Fit B-Splines with weighted penalization over differences of parameters
#'
#' @param XX_band The matrix \eqn{X^T X} where \code{X} is the design matrix. This argument is given
#' in the form of a band matrix, i.e., successive columns represent superdiagonals.
#' @param Xy The vector of currently estimated points \eqn{X^T y}, where \code{y} is the y-coordinate of the data.
#' @param degree The degree of the B-splines.
#' @param pen Positive penalty constant.
#' @param w Vector of weights. The case \eqn{\mathbf w = \mathbf 1} corresponds to fitting P-splines with difference #' order \code{degree + 1} (see \emph{Eilers, P., Marx, B. (1996) Flexible smoothing with B-splines and penalties}.)
#' @param old_par Initial parameter to serve as starting point of the iterating process.
#' @param maxiter Maximum number of Newton-Raphson iterations to be computed.
#' @param tol The tolerance chosen to diagnostic convergence of the adaptive ridge procedure.
#' @return The estimated parameter of the spline regression.
#' @export
wridge_solver <- function(XX_band, Xy, degree, pen,
                          w = rep(1, nrow(XX_band) - degree - 1),
                          old_par = rep(1, nrow(XX_band)),
                          maxiter = 1000,
                          tol = 1e-8) {
  for (iter in 1:maxiter) {
    par <- old_par + hessian_solver(old_par, XX_band, Xy,
                                    pen, w, diff = degree + 1)
    idx <- old_par != 0
    rel_error <- max(abs(par - old_par)[idx] / abs(old_par)[idx])
    if (rel_error < tol) break
    old_par <- par
  }
  par
}
#'
#' Fit B-splines with automatic knot selection.
#'
#' @param X design matrix
#' @param y Data y values
#' @param degree The degree of the splines. Recommended value is 3, which corresponds to natural splines.
#' @param pen A vector of positive penalty values. The adaptive spline regression is performed for every value of pen
#' @param maxiter Maximum number of iterations  in the main loop.
#' @param epsilon Value of the constant in the adaptive ridge procedure (see \emph{Frommlet, F., Nuel, G. (2016)
#' An Adaptive Ridge Procedure for L0 Regularization}.)
#' @param verbose Whether to print details at each step of the iterative procedure.
#' @param diff Order of the differences on the parameters. The value \code{degree + 1} is necessary to perform
#' selection of the knots.
#' @param tol The tolerance chosen to diagnostic convergence of the adaptive ridge procedure.
#' @export
aridge_solver <- function(X, y, pen, degree,
                          maxiter = 1000,
                          epsilon = 1e-5,
                          verbose = FALSE,
                          diff = degree + 1,
                          tol = 1e-6) {
  XX <- crossprod(X)
  XX_band <- cbind(bandsolve::mat2rot(XX + diag(rep(1e-20), ncol(X))), 0)
  Xy <- crossprod(X, y)
  # Define sigma0
  sigma0sq <- var(lm(y ~ X - 1)$residuals)
  # Define returned variables
  model <- X_sel <- knots_sel <- sel_ls <- par_ls <- vector("list", length(pen))
  aic <- bic <- ebic <- pen * NA
  dim <- loglik <- pen * NA
  # Initialize values
  old_sel <- rep(1, ncol(X) - diff)
  par <- rep(1, ncol(X))
  w <- rep(1, ncol(X) - diff)
  ind_pen <- 1
  # Main loop
  for (iter in 1:maxiter) {
    par <- wridge_solver(XX_band, Xy, degree,
                         pen[ind_pen], w,
                         old_par = par)
    w <- 1 / (diff(par, differences = diff) ^ 2 + epsilon ^ 2)
    sel <- w * diff(par, differences = diff) ^ 2
    if (verbose) {
      cat('iter =', iter, ' sum_sel = ', sum(sel), '\n')
      plot(sel, ylim = c(0, 1), main =
             cat('iter =', iter, ' sum_sel = ', sum(sel), '\n'))
    }
    converge <- max(abs(old_sel - sel)) < tol
    if (converge) {
      sel_ls[[ind_pen]] <- sel
      knots_sel[[ind_pen]] <- knots[sel > 0.99]
      X_sel[[ind_pen]] <- splines2::bSpline(
        x, knots = knots_sel[[ind_pen]], intercept = TRUE, degree = degree)
      model[[ind_pen]] <- lm(y ~ X_sel[[ind_pen]] - 1)
      if (verbose) {
        plot(x, y, col = "gray")
        lines(x, predict(model[[ind_pen]]), col = "red")
        lines(x, predict(lm(y ~ X - 1)), col = "blue")
        abline(v = knots_sel[[ind_pen]], col = "red")
      }
      par_ls[[ind_pen]] <- rep(NA, ncol(X))
      idx <- c(sel > 0.99, rep(TRUE, diff))
      par_ls[[ind_pen]][idx] <- model[[ind_pen]]$coefficients
      par_ls[[ind_pen]][!idx] <- 0
      loglik[ind_pen] <- log(sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq))
      dim[ind_pen] <- length(knots_sel[[ind_pen]]) + degree + 1
      bic[ind_pen] <- log(nrow(X)) * (length(knots_sel[[ind_pen]]) + degree + 1) +
        2 * log(sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq))
      aic[ind_pen] <- 2 * (length(knots_sel[[ind_pen]]) + degree + 1) +
        2 * log(sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq))
      ebic[ind_pen] <- bic[ind_pen] + 2 * lchoose(ncol(X), ncol(X_sel[[ind_pen]]))
      ind_pen <- ind_pen + 1
    }
    if (ind_pen > length(pen)) break
    old_sel <- sel
  }
  # Diagnostic for bad behavior of Adaptive Ridge
  sel_mat <- sel_ls %>%
    unlist() %>%
    round(digits = 1) %>%
    matrix(., ncol(X) - degree - 1)
  knots_sel_monotonous <- apply(sel_mat, 1, function(a) all(diff(a) <= 0))
  if (!all(knots_sel_monotonous)) {
    if (sum(!knots_sel_monotonous) >= 10) {
      warning(paste0("The models are not nested:\n",
                     sum(!knots_sel_monotonous),
                     " knots are dropped and then reselected"))
    } else {
    warning(paste0("The models are not nested:\n",
                   "Knots number ", paste(which(!knots_sel_monotonous), collapse = ', '),
                   " are dropped and then reselected"))
    }
  }
  # Print regularization path
  regul_df <- dplyr::data_frame(penalty = rep(pen, each = ncol(X)),
                         index = rep(1:(ncol(X)), length(pen)),
                         param = par_ls %>% unlist())
  path <- ggplot2::ggplot(regul_df, aes(penalty, param, color = as.factor(index))) +
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::geom_vline(xintercept = pen[which(diff(apply(sel_mat, 2, sum)) != 0) + 1],
               size = 0.2)
  # Return values
  list("sel" = sel_ls, "knots_sel" = knots_sel, "model" = model,
       "X_sel" = X_sel, "par" = par_ls, "sel_mat" = sel_mat,
       "aic" = aic, "bic" = bic, "ebic" = ebic, "path" = path,
       "dim" = dim, "loglik" = loglik)
}
#' @export
aspline <- function(x, y, knots = seq(min(x), max(x), length = 2 * length(x) + 2)[-c(1, 2 * length(x) + 2)],
                    pen = 10 ^ seq(-3, 6, length = 100),
                    degree = 3,
                    maxiter = 10000,
                    epsilon = 1e-5,
                    verbose = FALSE,
                    tol = 1e-6) {
  X <- splines2::bSpline(x, knots = knots, degree = degree, intercept = TRUE)
  XX <- crossprod(X)
  XX_band <- cbind(
    bandsolve::mat2rot(XX + diag(rep(1e-20), length(knots) + degree + 1)),
    0)
  Xy <- crossprod(X, y)
  # Define sigma0
  X_10 <- splines2::bSpline(
    x,
    knots = seq(min(x), max(x), length = 12)[-c(1, 12)],
    degree = degree,
    intercept = TRUE)
  sigma0sq <- var(lm(y ~ X_10 - 1)$residuals)
  rm(X)
  # Define returned variables
  model <- X_sel <- knots_sel <- sel_ls <- par_ls <- vector("list", length(pen))
  aic <- bic <- ebic <- loglik <- dim <- pen * NA
  # Initialize values
  old_sel <- rep(1, nrow(XX_band) - diff)
  par <- rep(1, nrow(XX_band))
  w <- rep(1, nrow(XX_band) - degree - 1)
  ind_pen <- 1
  # Main loop
  for (iter in 1:maxiter) {
    par <- wridge_solver(XX_band, Xy, degree,
                         pen[ind_pen], w,
                         old_par = par)
    w <- 1 / (diff(par, differences = diff) ^ 2 + epsilon ^ 2)
    sel <- w * diff(par, differences = diff) ^ 2
    if (verbose) {
      cat('iter =', iter, ' sum_sel = ', sum(sel), '\n')
      plot(sel, ylim = c(0, 1), main =
             cat('iter =', iter, ' sum_sel = ', sum(sel), '\n'))
    }
    converge <- max(abs(old_sel - sel)) < tol
    if (converge) {
      sel_ls[[ind_pen]] <- sel
      knots_sel[[ind_pen]] <- knots[sel > 0.99]
      X_sel[[ind_pen]] <- splines2::bSpline(
        x, knots = knots_sel[[ind_pen]], intercept = TRUE, degree = degree)
      model[[ind_pen]] <- lm(y ~ X_sel[[ind_pen]] - 1)
      # if (verbose) {
      #   plot(x, y, col = "gray")
      #   lines(x, predict(model[[ind_pen]]), col = "red")
      #   lines(x, predict(lm(y ~ X - 1)), col = "blue")
      #   abline(v = knots_sel[[ind_pen]], col = "red")
      # }
      par_ls[[ind_pen]] <- rep(NA, nrow(XX_band))
      idx <- c(sel > 0.99, rep(TRUE, diff))
      par_ls[[ind_pen]][idx] <- model[[ind_pen]]$coefficients
      par_ls[[ind_pen]][!idx] <- 0
      loglik[ind_pen] <- log(sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq))
      dim[ind_pen] <- length(knots_sel[[ind_pen]]) + degree + 1
      bic[ind_pen] <- log(ncol(XX_band)) * (length(knots_sel[[ind_pen]]) + degree + 1) +
        2 * log(sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq))
      aic[ind_pen] <- 2 * (length(knots_sel[[ind_pen]]) + degree + 1) +
        2 * log(sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq))
      ebic[ind_pen] <- bic[ind_pen] + 2 * lchoose(length(knots) + degree + 1, ncol(X_sel[[ind_pen]]))
      ind_pen <- ind_pen + 1
    }
    if (ind_pen > length(pen)) break
    old_sel <- sel
  }
  # Diagnostic for bad behavior of Adaptive Ridge
  sel_mat <- sel_ls %>%
    unlist() %>%
    round(digits = 1) %>%
    matrix(length(knots))
  knots_sel_monotonous <- apply(sel_mat, 1, function(a) all(diff(a) <= 0))
  if (!all(knots_sel_monotonous)) {
    if (sum(!knots_sel_monotonous) >= 10) {
      warning(paste0("The models are not nested: ",
                     sum(!knots_sel_monotonous),
                     " knots are dropped and then reselected"))
    } else {
    warning(paste0("The models are not nested: ",
                   "knots number ", paste(which(!knots_sel_monotonous), collapse = ', '),
                   " are dropped and then reselected"))
    }
  }
  # Print regularization path
  regul_df <- dplyr::data_frame(penalty = rep(pen, each = nrow(XX_band)),
                                index = rep(1:(nrow(XX_band)), length(pen)),
                                param = par_ls %>% unlist())
  path <- ggplot2::ggplot(regul_df, aes(penalty, param, color = as.factor(index))) +
    ggplot2::geom_step() +
    ggplot2::scale_x_log10() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::geom_vline(xintercept = pen[which(diff(apply(sel_mat, 2, sum)) != 0) + 1],
                        size = 0.2)
  # Return values
  list("model" = model, "sel" = sel_ls, "knots" = knots_sel,
       "design" = X_sel, "par" = par_ls, "sel_mat" = sel_mat,
       "aic" = aic, "bic" = bic, "ebic" = ebic, "loglik" = loglik,
       "dim" = dim, "path" = path)
}
#' K-fold cross-validation
#' @export
kcv <- function(x, y, pen = 10 ^ seq(-3, 3, length = 50), nfold = 10) {
  x <- as.vector(x)
  y <- as.vector(y)
  score_matrix <- matrix(NA, nfold, length(pen))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(length(x)/nfold) * (ind - 1) + 1,
                       floor(length(x)/nfold) * ind)
    sample_train <- setdiff(1:length(x), sample_test)
    x_train <- x[sample_train]
    y_train <- y[sample_train]
    x_test <- x[sample_test]
    y_test <- y[sample_test]
    train <- aspline(x_train, y_train, pen)$model
    score_matrix[ind, ] <- sapply(
      seq_along(train),
    function(model) log(sum((y_test - predict(model, x_test)) ^ 2))
    )
    if (any(is.null(score_matrix[ind, ]))) {
      stop('Error in call to aspline')
    }
  }
  colSums(score_matrix)
}
