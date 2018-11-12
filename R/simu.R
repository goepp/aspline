#' @export
hetero_1 <- function(x, alpha_var = 1) {
  (x * 0.3 + sqrt(alpha_var) * 0.2) ^ 2 * alpha_var
}
#' @export
logit <- function(x) {
  1 / (1 + exp(-20 * (x - 0.5)))
}
#' @export
bump <- function(x) {
  (x + 2 * exp(-(16 * (x - 0.5)) ^ 2)) / 2.5
}
#' @export
spa_het <- function(x) {
  sqrt(x * (1 - x)) * sin((2 * pi * (1 + 2 ^ (-3 / 5))) / (x + 2 ^ (-3 / 5))) + 0.5
}
#' @export
sine3 <- function(x) {
  sin(3 * 2 * pi * x) * 0.5 + 0.5
}
#' @export
sine12 <- function(x) {
  sin(12 * 2 * pi * x) / 2 + 0.5
}
#' @export
compute_mse <- function(x, y, k, method) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    s_spline <- pspline::sm.spline(x, y)
    # lines(x, predict(s_spline, x))
    sum((y - predict(s_spline, x)) ^ 2) / length(y)
  } else if (method == "p") {
    p_spline <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
    # lines(x, predict(p_spline), col = "green")
    sum(p_spline$residuals ^ 2) / length(y)
  } else if (method %in% c("a", "a_aic")) {
    a_spline <- aridge_solver(x, y, knots)
    sum(a_spline$fit$aic$residuals ^ 2) / length(y)
  } else if (method == "a_bic") {
    a_spline <- aridge_solver(x, y, knots)
    sum(a_spline$fit$bic$residuals ^ 2) / length(y)
  } else if (method == "a_ebic") {
    a_spline <- aridge_solver(x, y, knots)
    sum(a_spline$fit$ebic$residuals ^ 2) / length(y)
  } else if (method == "bars") {
    bars <- barsN.fun(x, y, priorparam = c(1, 20))
    sum((bars$postmodes - y) ^ 2 / length(y))
  } else if (method == "fks") {
    fks <- fit.search.numknots(x, y, degree = 3, search = "genetic")
    sum((fitted.freekt(fks) - y) ^ 2 / length(y))
  } else if (method == "freeps") {
    freeps <- freepsgen(x, y, degree = 3, numknot = 5)
    sum((fitted.freekt(freeps) - y) ^ 2 / length(y))
  } else {
    stop("Error: method argument not correct")
  }
}
#' @export
fitted_wrapper <- function(x, y, k, method) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    s_spline <- pspline::sm.spline(x, y)
    as.vector(s_spline$ysmth)
  } else if (method == "p") {
    p_spline <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
    p_spline$fitted.values
  } else if (method %in% c("a", "a_aic")) {
    a_spline <- aridge_solver(x, y, knots)
    a_spline$fit$aic$fitted.values
  } else if (method == "a_bic") {
    a_spline <- aridge_solver(x, y, knots)
    a_spline$fit$bic$fitted.values
  } else if (method == "a_ebic") {
    a_spline <- aridge_solver(x, y, knots)
    a_spline$fit$ebic$fitted.values
  } else if (method == "bars") {
    bars <- barsN.fun(x, y, priorparam = c(1, 20))
    bars$postmodes
  } else if (method == "fks") {
    fks <- fit.search.numknots(x, y, degree = 3, search = "genetic")
    as.vector(fitted.freekt(fks))
  } else if (method == "freeps") {
    freeps <- freepsgen(x, y, degree = 3, numknot = 5)
    as.vector(fitted.freekt(freeps))
  } else {
    stop("Error: method argument not correct")
  }
}
#' @export
difference_s <- function(x, fit, fun) {
  (as.vector(predict(fit, x = x)) - pryr::fget(fun)(x)) ^ 2
}
#' @export
difference_p <- function(x, fit, fun) {
  (as.vector(predict(fit, data_frame(x = x))) - pryr::fget(fun)(x)) ^ 2
}
#' @export
difference_a <- function(x, fit, fun) {
  (as.vector(predict(fit, data_frame(x = x))) - pryr::fget(fun)(x)) ^ 2
}
#' @export
l2_fit <- function(x, y, k, method, fun) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    fit <- pspline::sm.spline(x, y)
    error <- integrate(difference_s, 0, 1, fit = fit, fun = fun)$value
  } else if (method == "p") {
    fit <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
    error <- integrate(difference_p, 0, 1, fit = fit, fun = fun)$value
  } else if (method %in% c("a", "a_aic")) {
    aridge <- aridge_solver(x, y, knots)
    fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$aic)]]))
    error <- integrate(difference_a, 0, 1, fit = fit, fun = fun)$value
  } else if (method == "a_bic") {
    aridge <- aridge_solver(x, y, knots)
    fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$bic)]]))
    error <- integrate(difference_a, 0, 1, fit = fit, fun = fun)$value
  } else if (method == "a_ebic") {
    aridge <- aridge_solver(x, y, knots)
    fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$ebic)]]))
    error <- integrate(difference_a, 0, 1, fit = fit, fun = fun)$value
  } else if (method == "bars") {
    bars <- barsN.fun(x, y, priorparam = c(1, 20))
    fit <- approxfun(x, bars$postmodes)
    error <- integrate(function(x) (fit(x) - pryr::fget(fun)(x)) ^ 2, 0, 1)$value
  } else if (method == "fks") {
    fks <- fit.search.numknots(x, y, degree = 3, search = "genetic")
    fit <- approxfun(x, fitted.freekt(fks))
    error <- integrate(function(x) (fit(x) - pryr::fget(fun)(x)) ^ 2, 0, 1)$value
  } else if (method == "freeps") {
    freeps <- freepsgen(x, y, degree = 3, numknot = 5)
    fit <- approxfun(x, fitted.freekt(freeps))
    error <- integrate(function(x) (fit(x) - pryr::fget(fun)(x)) ^ 2, 0, 1)$value
  } else {
    stop("Error: method argument not correct")
  }
  error
}
#' @export
predict_wrapper <- function(x, y, k, x_seq, method) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    fit <- pspline::sm.spline(x, y)
    as.vector(predict(fit, x = x_seq))
  } else if (method == "p") {
    fit <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
    as.vector(predict(fit, data_frame(x = x_seq)))
  } else if (method %in% c("a", "a_aic")) {
    aridge <- aridge_solver(x, y, knots)
    fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$aic)]]))
    as.vector(predict(fit, data_frame(x = x_seq)))
  } else if (method == "a_bic") {
    aridge <- aridge_solver(x, y, knots)
    fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$bic)]]))
    as.vector(predict(fit, data_frame(x = x_seq)))
  } else if (method == "a_ebic") {
    aridge <- aridge_solver(x, y, knots)
    fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$ebic)]]))
    as.vector(predict(fit, data_frame(x = x_seq)))
  } else {
    stop("Error: method argument not correct")
  }
}
#' @export
boldify <- function(mat, ind_row = 1:ncol(mat)) {
  ind_bool <-  1:ncol(mat) %in% ind_row
  submat <- mat[, ind_row]
  for (row in 1:nrow(submat)) {
    ind_max <- which.min(submat[row, ])
    submat[row, ind_max] <- paste0("$\\bm{", submat[row, ind_max], "}$")
  }
  result <- matrix("", nrow(mat), ncol(mat))
  result[, which(ind_bool)] <- as.matrix(submat)
  result[, !ind_bool] <- mat[, !ind_bool]
  result %>%
    "colnames<-"(colnames(mat)) %>%
    "rownames<-"(rep("", nrow(mat)))
}
#' @export
nlevel_fit <- function(x, y, k, method, fun) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    nlevel <- NA
  } else if (method == "p") {
    nlevel <- NA
  } else if (method %in% c("a", "a_aic")) {
    aridge <- aridge_solver(x, y, knots)
    # fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$aic)]]))
    nlevel <- ncol(bs(x, knots = aridge$knots_sel[[which.min(aridge$aic)]]))
  } else if (method == "a_bic") {
    aridge <- aridge_solver(x, y, knots)
    # fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$bic)]]))
    nlevel <- ncol(bs(x, knots = aridge$knots_sel[[which.min(aridge$bic)]]))
  } else if (method == "a_ebic") {
    aridge <- aridge_solver(x, y, knots)
    # fit <- lm(y ~ bs(x, knots = aridge$knots_sel[[which.min(aridge$ebic)]]))
    nlevel <- ncol(bs(x, knots = aridge$knots_sel[[which.min(aridge$ebic)]]))
  } else if (method == "bars") {
    bars <- barsN.fun(x, y, priorparam = c(1, 20))
    nlevel <- mode(bars$no.knots)
    } else {
    stop("Error: method argument not correct")
  }
  nlevel
}
