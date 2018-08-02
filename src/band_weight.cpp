#include <Rcpp.h>
using namespace Rcpp;

//' Create the penalty matrix
//'
//' @param w Vector of weights
//' @param diff Order of the differences to be applied to the parameters. Must be a strictly positive integer
//' @return Weighted penalty matrix \eqn{D^T diag(w) D} where
//'  \code{D <- diff(diag(length(w) + diff), differences = diff)}. Only the non-null superdiagonals of
//'  the weight matrix are returned, each column corresponding to a diagonal.
//' @export
// [[Rcpp::export]]
NumericMatrix band_weight(NumericVector w, int diff) {
  int ws = w.size();
  // Compute the entries of the difference matrix
  NumericVector binom(diff + 1);
  for (int i = 0; i <= diff; i ++) {
    binom(i) = Rf_choose(diff, i) * pow((-1), i);
  }
  // Compute the limit indices
  NumericMatrix ind_mat(ws + diff, 2);
  for (int ind = 0; ind < ws + diff; ind ++) {
    ind_mat(ind, 0) = ind - diff < 0 ? 0 : ind - diff;
    ind_mat(ind, 1) = ind < ws - 1 ? ind : ws - 1;
  }
  // Main loop
  NumericMatrix result(ws + diff, diff + 1);
  for (int j = 0; j < result.ncol(); j ++) {
    for (int i = 0; i < ws + diff - j; i ++) {
      double temp = 0.;
      for (int k = ind_mat(i + j, 0); k <= ind_mat(i, 1); k ++) {
        temp += binom(i - k) * binom(i + j - k) * w(k);
      }
      result(i, j) = temp;
    }
  }
  return result;
}


