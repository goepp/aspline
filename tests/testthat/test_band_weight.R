context("band_weight.R")

library(Rcpp)

diff <- 10
w <- rnorm(1000)

band_weight(w, diff)
