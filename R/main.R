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
  bandsolve(XX_band + pen * band_weight(w, diff), Xy) - par
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
#' @param x,y Input data, numeric vectors of same length
#' @param knots Knots
#' @param pen A vector of positive penalty values. The adaptive spline regression is performed for every value of pen
#' @param degree The degree of the splines. Recommended value is 3, which corresponds to natural splines.
#' @param family A description of the error distribution and link function to be used in the model. The "gaussian", "binomial", and "poisson" families are currently implemented, corresponding to the linear regression, logistic regression, and Poisson regression, respectively.
#' @param maxiter Maximum number of iterations  in the main loop.
#' @param epsilon Value of the constant in the adaptive ridge procedure (see \emph{Frommlet, F., Nuel, G. (2016)
#' An Adaptive Ridge Procedure for L0 Regularization}.)
#' @param verbose Whether to print details at each step of the iterative procedure.
#' @param tol The tolerance chosen to diagnostic convergence of the adaptive ridge procedure.
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{sel}: list giving for each value of \code{lambda} the vector of the knot selection weights (a knot is selected if its weight is equal to 1.)}
#' \item{\code{knots_sel}: list giving for each value of \code{lambda} the vector of selected knots.}
#' \item{\code{model}: list giving for each value of \code{lambda} the fitted regression model.}
#' \item{\code{par}: parameters of the models for each value of \code{lambda}.}
#' \item{\code{sel_mat}: matrix of booleans whose columns indicate whether each knot is selected.}
#' \item{\code{aic}, \code{bic}, and \code{ebic}: Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), and Extended BIC (EBIC) scores, for each value of \code{lambda}.}
#' \item{\code{dim}: number of selected knots for each value of \code{lambda}.}
#' \item{\code{loglik}: log-likelihood of the selected model, for each value of \code{lambda}.}
#' }
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom stats var
#' @importFrom stats lm
#' @importFrom stats predict
#' @importFrom rlang .data
#' @export
aspline <- function(x, y,
                    knots = seq(min(x), max(x), length = 42)[-c(1, 42)],
                    pen = 10 ^ seq(-3, 3, length = 100),
                    degree = 3L,
                    family = c("gaussian", "binomial", "poisson"),
                    maxiter = 1000,
                    epsilon = 1e-5,
                    verbose = FALSE,
                    tol = 1e-6) {
  family <- match.arg(family)
  X <- splines2::bSpline(x, knots = knots, intercept = TRUE, degree = degree)
  if (family == "gaussian") {
    XX <- crossprod(X)
    XX_band <- cbind(mat2rot(XX + diag(rep(1e-20), ncol(X))), 0)
    Xy <- crossprod(X, y)
  } else  if (family %in% c("binomial", "poisson")) {
    comp <- block_design(X, degree)
    B <- comp$B
    alpha <- comp$alpha
  }
  # Define sigma0
  sigma0sq <- var(lm(y ~ X - 1)$residuals)
  # Define returned variables
  model <- X_sel <- knots_sel <- sel_ls <- par_ls <- vector("list", length(pen))
  aic <- bic <- ebic <- pen * NA
  dim <- loglik <- pen * NA
  # Initialize values
  old_sel <- rep(1, ncol(X) - degree - 1)
  par <- rep(1, ncol(X))
  w <- rep(1, ncol(X) - degree - 1)
  ind_pen <- 1

  # Main loop
  for (iter in 1:maxiter) {
    if (family == "gaussian") {
      par <- wridge_solver(XX_band, Xy, degree,
                           pen[ind_pen], w,
                           old_par = par)
    } else if (family %in% c("binomial", "poisson")) {
      par <- wridge_solver_glm(X, y, B, alpha, degree,pen[ind_pen],
                               family = family,
                               old_par = par,
                               w = w,
                               maxiter = 1000)$par
    }
    w <- 1 / (diff(par, differences = degree + 1) ^ 2 + epsilon ^ 2)
    sel <- w * diff(par, differences = degree + 1) ^ 2
    converge <- max(abs(old_sel - sel)) < tol
    if (converge) {
      sel_ls[[ind_pen]] <- sel
      knots_sel[[ind_pen]] <- knots[sel > 0.99]
      design <- splines2::bSpline(
        x, knots = knots_sel[[ind_pen]], intercept = TRUE, degree = degree)
      X_sel[[ind_pen]] <- design
      model[[ind_pen]] <- lm(y ~ design - 1)
      if (verbose) {
        plot(x, y, col = "gray")
        lines(x, predict(model[[ind_pen]]), col = "red")
        lines(x, predict(lm(y ~ X - 1)), col = "blue")
        abline(v = knots_sel[[ind_pen]], col = "red")
      }
      par_ls[[ind_pen]] <- rep(NA, ncol(X))
      idx <- c(sel > 0.99, rep(TRUE, degree + 1))
      par_ls[[ind_pen]][idx] <- model[[ind_pen]]$coefficients
      par_ls[[ind_pen]][!idx] <- 0
      loglik[ind_pen] <- 1 / 2 * sum((model[[ind_pen]]$residuals) ^ 2 / sigma0sq)
      dim[ind_pen] <- length(knots_sel[[ind_pen]]) + degree + 1
      aic[ind_pen] <- 2 * dim[ind_pen] + 2 * loglik[ind_pen]
      bic[ind_pen] <- log(nrow(X)) * dim[ind_pen] + 2 * loglik[ind_pen]
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
    matrix(ncol(X) - degree - 1)
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
  # Return values
  list("sel" = sel_ls, "knots_sel" = knots_sel, "model" = model,
       "par" = par_ls, "sel_mat" = sel_mat,
       "aic" = aic, "bic" = bic, "ebic" = ebic,
       "dim" = dim, "loglik" = loglik)
}
#' @export
#' @describeIn aspline Alias for \code{aspline}, for backwards compatibility.
aridge_solver <- aspline
hessian_solver_glm <- function(par, X, y, degree, pen, family,
                               w = rep(1, ncol(X) - degree - 1)) {
  if (family == "gaussian") {
    W <- diag(length(par))
    g_inv <- identity
    g_p <- function(x) return(0)
  }
  if (family == "poisson") {
    g_inv <- exp
    g_p <- function(x) 1 / x
    W <- diag(g_inv(as.vector(X %*% par)))
  }
  if (family == "binomial") {
    temp <- as.vector(exp(-X %*% par))
    W <- diag(temp / (1 + temp) ^ 2)
    g_inv <- function(x) 1 / (1 + exp(-x))
    g_p <- function(x) 1 / ((1 - x) * x)
  }
  D <- diff(diag(ncol(X)), differences = degree + 1)
  mat <- t(X) %*% W %*% X  + pen * t(D) %*% diag(w) %*% D
  # vect <- t(X) %*% W %*% (y - X %*% par + g_inv(X %*% par)) # OLD AND FALSE
  vect <- t(X) %*% (W %*% X %*% par + y - g_inv(X %*% par))
  as.vector(solve(mat, vect))
}
hessian_solver_glm_band <- function(par, X, y, B, alpha, pen, w, degree,
                                    family = c("gaussian", "binomial", "poisson")) {
  family <- match.arg(family)
  if (family == "gaussian") {
    glm_weight <- rep(1, length(par))
    g_inv <- identity
    g_p <- function(x) return(0)
  }
  if (family == "poisson") {
    glm_weight <- g_inv(as.vector(X %*% par))
    g_inv <- exp
    g_p <- function(x) 1 / x
  }
  if (family == "binomial") {
    temp <- as.vector(exp(-X %*% par))
    glm_weight <- temp / (1 + temp) ^ 2
    g_inv <- function(x) 1 / (1 + exp(-x))
    g_p <- function(x) 1 / ((1 - x) * x)
  }
  # if (family == "gaussian") glm_weight <- rep(1, length(par))
  # if (family == "poisson") glm_weight <- as.vector(X %*% par)
  # if (family == "binomial") {
  #   temp <- as.vector(exp(-X %*% par))
  #   glm_weight <- temp / (1 + temp) ^ 2
  # }
  XWX_band <- cbind(weight_design_band(glm_weight, alpha, B), 0)
  mat <- XWX_band + pen * band_weight(w, degree + 1)
  # vect <- sweep(t(X), MARGIN = 2, glm_weight, `*`) %*% y
  vect <- crossprod(X, sweep(X, 1, glm_weight, `*`) %*% par + y - g_inv(X %*% par))
  as.vector(bandsolve(mat, vect))
}
wridge_solver_glm <- function(X, y, B, alpha, degree, pen,
                              family = c("gaussian", "poisson", "binomial"),
                              old_par = rep(1, ncol(X)),
                              w = rep(1, ncol(X) - degree - 1),
                              maxiter = 1000) {
  family <- match.arg(family)
  for (iter in 1:maxiter) {
    # as.vector((X %*% old_par)) %>% plot()
    # par <- hessian_solver_glm_band(old_par, X, y, B, alpha, pen, w, degree, family = family)
    par <- hessian_solver_glm(old_par, X, y, degree, pen, family, w)
    # (par) %>% plot()
    idx <- old_par != 0
    rel_error <- max(abs(par - old_par)[idx] / abs(old_par)[idx])
    if (rel_error < 1e-5) break
    old_par <- par
  }
  if (iter == maxiter) warnings("Warning: NR did not converge.")
  list('par' = par, 'iter' = iter)
}


