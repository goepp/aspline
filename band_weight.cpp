#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' Create the matrix of weights
//'
//' @param w The vector of weights
//' @param diff The order of the differences to be applied to the parameters. Must be a strictly positive integer
//' @return The matrix of weights \code{t(D) %*% diag(w) %*% D} where
//'  \code{D <- diff(diag(length(w) + diff), differences = diff)}. Only the non-null superdiagonals of
//'  the weight matrixare returned, each column corresponding to a diagonal.
//' @export
// [[Rcpp::export]]
// [[Rcpp::interfaces(r, cpp)]]
NumericMatrix band_weight(NumericVector w, int diff) {
  int ws = w.size();
  NumericVector binom(diff + 1);
  for (int i = 0; i <= diff; i ++) {
    binom(i) = Rf_choose(diff, i) * pow((-1), i);
  }
  NumericMatrix ind_mat(ws + diff, 2);
  for (int ind = 0; ind < ws + diff; ind ++) {
    ind_mat(ind, 0) = ind - diff < 0 ? 0 : ind - diff;
    ind_mat(ind, 1) = ind < ws - 1 ? ind : ws - 1;
  }
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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
band_weight(rnorm(10), 2)

*/
