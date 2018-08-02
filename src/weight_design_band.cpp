#include <Rcpp.h>
using namespace Rcpp;

//' Fast computation of weighted design matrix for generalized linear model
//'
//' @param w Vector of weights.
//' @param alpha Vector of indexes representing the start of blocks of the design matrix, as given by \link{block_design}.
//' @param B Design matrix in compressed block format, as given by \link{block_design}.
//' @return Weighted design matrix \eqn{X^T diag(w) X} where \code{X} is the design matrix and \code{W = diag(w)} is
//' a diagonal matrix of weights.
//' @export
// [[Rcpp::export]]
NumericMatrix weight_design_band(NumericVector w, NumericVector alpha, NumericMatrix B) {
  int as = alpha.size();
  int m = B.ncol();
  NumericMatrix ind_mat(as + m - 2, 2);
  for (int ind = 0; ind < as + m - 2; ind ++) {
    ind_mat(ind, 0) = ind - m + 1 < 0 ? 0 : ind - m + 1;
    ind_mat(ind, 1) = ind < as - m ? ind : as - m;
  }
  NumericMatrix M(as + m - 2, m);
  for (int k = 0; k < m; k ++) {
    for (int i = 0; i < as + m - 2 - k; i ++) {
      int u_min = ind_mat(i, 0) < ind_mat(k + i, 0) ? ind_mat(k + i, 0) : ind_mat(i, 0);
      int u_max = ind_mat(i, 1) < ind_mat(k + i, 1) ? ind_mat(i, 1) : ind_mat(k + i, 1);
      double result = 0.;
      for (int u = u_min; u <= u_max; u ++) {
        int l_max = alpha(u + 1) - alpha(u) - 1;
        for (int l = 0; l <= l_max; l ++) {
          result += B(alpha(u) + l - 1, i - u) * B(alpha(u) + l - 1, i + k - u) * w(alpha(u) + l - 1);
        }
      }
      M(i, k) = result;
    }
  }
  return M;
}

