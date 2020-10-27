hetero_1 <- function(x, alpha_var = 1) {
  (x * 0.3 + sqrt(alpha_var) * 0.2) ^ 2 * alpha_var
}
logit <- function(x) {
  1 / (1 + exp(-20 * (x - 0.5)))
}
bump <- function(x) {
  (x + 2 * exp(-(16 * (x - 0.5)) ^ 2)) / 2.5
}
doppler <- function(x) {
  ifelse(x > 0, sin(4 / x ^ 0.7) + 1.5, 0.5)
}
wave <- function(x) {
  (sin(4 * pi * x ^ 4) * (1 - x) + 0.22) / 0.55
}
spa_het <- function(x) {
  sqrt(x * (1 - x)) * sin((2 * pi * (1 + 2 ^ (-3 / 5))) / (x + 2 ^ (-3 / 5))) + 0.5
}
sine3 <- function(x) {
  sin(3 * 2 * pi * x) * 0.5 + 0.5
}
sine12 <- function(x) {
  sin(12 * 2 * pi * x) / 2 + 0.5
}
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
    # dyn.load("barsN.so", now = F)
    bars <- barsN.fun(x, y, priorparam = c(1, k))
    bars$postmodes
  } else if (method == "fks") {
    fks <- freeknotsplines::fit.search.numknots(x, y, degree = 3,
                                                maxknot = k, search = "genetic")
    as.vector(fitted.freekt(fks))
  } else if (method == "freeps") {
    freeps <- freeknotsplines::freepsgen(x, y, degree = 3, numknot = 5)
    as.vector(fitted.freekt(freeps))
  } else {
    stop("Error: method argument not correct")
  }
}
difference_s <- function(x, fit, fun) {
  (as.vector(predict(fit, x = x)) - pryr::fget(fun)(x)) ^ 2
}
difference_p <- function(x, fit, fun) {
  (as.vector(predict(fit, data_frame(x = x))) - pryr::fget(fun)(x)) ^ 2
}
difference_a <- function(x, fit, fun) {
  (as.vector(predict(fit, data_frame(x = x))) - pryr::fget(fun)(x)) ^ 2
}
difference_a_new <- function(x, param, knots, fun) {
  B <- fda::bsplineS(x = x, breaks = knots, norder = 3 + 1, returnMatrix = TRUE)
  (as.vector(B %*% param - pryr::fget(fun)(x))) ^ 2
}
error_wrapper <- function(ind, design, data) {
  sample <- data %>% filter(ind_wrapper == ind)
  l2_fit(sample$x, sample$y, sample$k[1], sample$method[1], sample$fun[1])
}
#' @importFrom stats integrate
#' @importFrom stats rnorm
#' @importFrom stats approxfun
l2_fit <- function(x, y, k, method, fun) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    fit <- pspline::sm.spline(x, y)
    error <- integrate(difference_s, 0, 1, fit = fit, fun = fun)$value
  } else if (method == "p") {
    fit <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
    error <- integrate(difference_p, 0, 1, fit = fit, fun = fun)$value
  } else if (method == "a_new") {
    temp <- aspline(x, y, degree = 3, n_knots = k)
    param <- temp$result[which.min(temp$cv), ]
    knots <- seq(min(x) - diff(range(x)) / (k - 6),
                 max(x) + diff(range(x))  / (k - 6), length.out = k)
    error <- integrate(difference_a_new, 0, 1, param = param, knots = knots, fun = fun)$value
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
    dyn.load("barsN.so", now = F)
    bars <- barsN.fun(x, y, priorparam = c(1, length(knots)))
    fit <- approxfun(x, bars$postmodes)
    error <- integrate(function(x) (fit(x) - pryr::fget(fun)(x)) ^ 2, 0, 1)$value
  } else if (method == "fks") {
    fks <- freeknotsplines::fit.search.numknots(x, y, degree = 3, search = "genetic",
                                                maxknot = length(knots))
    fit <- approxfun(x, fitted.freekt(fks))
    error <- integrate(function(x) (fit(x) - pryr::fget(fun)(x)) ^ 2, 0, 1)$value
  } else if (method == "freeps") {
    freeps <- freeknotsplines::freepsgen(x, y, degree = 3, numknot = 5)
    fit <- approxfun(x, fitted.freekt(freeps))
    error <- integrate(function(x) (fit(x) - pryr::fget(fun)(x)) ^ 2, 0, 1)$value
  } else {
    stop("Error: method argument not correct")
  }
  error
}
nlevel_wrapper <- function(ind, design, data) {
  sample <- data %>% dplyr::filter("ind_wrapper" == ind)
  c(ind_wrapper = ind,
    nlevel = nlevel_fit(sample$x, sample$y, sample$k[1], sample$method[1], sample$fun[1]))
}
predict_wrapper <- function(x, y, k, x_seq, method) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    fit <- pspline::sm.spline(x, y)
    as.vector(predict(fit, x = x_seq))
  } else if (method == "p") {
    fit <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
    as.vector(predict(fit, data_frame(x = x_seq)))
  } else if (method == "a_new") {
    temp <- aspline(x, y, degree = 3, n_knots = k)
    param <- temp$result[which.min(temp$cv), ]
    knots <- seq(min(x) - diff(range(x)) / (k - 6),
                 max(x) + diff(range(x))  / (k - 6), length.out = k)
    B <- fda::bsplineS(x = x_seq, breaks = knots, norder = 3 + 1, returnMatrix = TRUE)
    as.vector(B %*% param)
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
nknot_fit <- function(x, y, k, method, degree) {
  knots <- seq(0, 1, length = k + 2)[-c(1, k + 2)]
  if (method == "s") {
    nknot <- NA
  } else if (method == "p") {
    nknot <- NA
  } else if (method %in% c("a", "a_aic")) {
    aridge <- aridge_solver(x, y, knots = knots, degree = degree)
    nknot <- length(aridge$knots_sel[[which.min(aridge$aic)]])
  } else if (method == "a_bic") {
    aridge <- aridge_solver(x, y, knots = knots, degree = degree)
    nknot <- length(aridge$knots_sel[[which.min(aridge$bic)]])
  } else if (method == "a_ebic") {
    aridge <- aridge_solver(x, y, knots = knots, degree = degree)
    nknot <- length(aridge$knots_sel[[which.min(aridge$ebic)]])
  } else if (method == "bars") {
    if (degree != 3) warning("BARS is only available with degree = 3")
    bars <- barsN.fun(x, y, priorparam = c(1, length(knots)))
    ux <- unique(bars$no.knots)
    nknot <- ux[which.max(tabulate(match(bars$no.knots, ux)))]
  } else if (method == "fks") {
    fks <- freeknotsplines::fit.search.numknots(x, y, degree = degree,
                                                search = "genetic")
    nknot <- length(fks@optknot)
  } else {
    stop("Error: method argument not correct")
  }
  return(nknot)
}
gen_data_hetero <- function(ind, design) {
  sample <- data_frame(
    x = seq(0, 1, length = design$sample_size[ind]),
    y = pryr::fget(design$fun[ind])(x) + rnorm(length(x), 0, sd = hetero_1(x, design$sigma[ind]))
  ) %>%
    mutate(sample_size = design$sample_size[ind]) %>%
    mutate(sigma = design$sigma[ind]) %>%
    mutate(ind_rep = design$ind_rep[ind]) %>%
    mutate(k = design$k[ind]) %>%
    mutate(fun = design$fun[ind]) %>%
    mutate(method = design$method[ind]) %>%
    mutate(ind_wrapper = ind)
  return(sample)
}
gen_data <- function(ind, design) {
  sample <- data_frame(
    x = seq(0, 1, length = design$sample_size[ind]),
    y = pryr::fget(design$fun[ind])(x) + rnorm(length(x), 0, sd = design$sigma[ind])
  ) %>%
    mutate(sample_size = design$sample_size[ind]) %>%
    mutate(sigma = design$sigma[ind]) %>%
    mutate(ind_rep = design$ind_rep[ind]) %>%
    mutate(k = design$k[ind]) %>%
    mutate(fun = design$fun[ind]) %>%
    mutate(method = design$method[ind]) %>%
    mutate(degree = design$degree[1]) %>%
    mutate(ind_wrapper = ind)
  sample
}
