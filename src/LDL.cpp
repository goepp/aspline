#include <Rcpp.h>
using namespace Rcpp;

//' @description Fast inplace LDL decomposition of symmetric band matrix of length k.
//' @title LDL
//' @param D Rotated row-wised matrix of dimensions n*k, with first column corresponding to the diagonal, the second to the first super-diagonal and so on.
//' @return List with D as solution of our LDL decomposition.
//' @name LDL
//' @examples 
//' 
//' n=10;
//' D0=1:10;
//' D1=exp(-c(1:9));
//' D=cbind(D0,c(D1,0))
//' sol=LDL(D)
//' @export

// [[Rcpp::export]]
List LDL(NumericMatrix D) {
  int n = D.nrow();
  int K = D.ncol()-1;
  
// LDL in-place 
  for (int i=1; i<=n; i++) {
    int j0=i-K;
    if (j0<1) j0=1;
    for (int j=j0; j<=i; j++) {
      for (int k=j0; k<j; k++)
      D(j-1,i-j)-=D(k-1,i-k)*D(k-1,j-k)*D(k-1,0);
      if (i>j) D(j-1,i-j)/=D(j-1,0);
    }
  }
  return Rcpp::List::create(Rcpp::Named("D")=D);
}