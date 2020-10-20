#include <Rcpp.h>
#include "LDL.h"
using namespace Rcpp;

// [[Rcpp::export]]
List bandsolve_cpp(NumericMatrix D, NumericMatrix B) {
    int n = D.nrow();
    LDL(D);
    int K = D.ncol()-1;
    
    int x1 = B.ncol();
    
    for (int l=0; l<x1;l++) { //solve over each column
      
    for (int i=2; i<=n; i++) { //solve b=inv(L)b
    int jmax=i-1;
    if (jmax>K) jmax=K;
    for (int j=1; j<=jmax; j++)
    B(i-1,l)-=D(i-j-1,j)*B(i-j-1,l);
  }
    // solve b=b/D
  for (int i=0; i<n; i++) B(i,l)/=D(i,0);
    // solve b=inv(t(L))b=inv(L*D*t(L))b
  for (int i=n-1; i>=1; i--) {
    int jmax=n-i;
    if (jmax>K) jmax=K;
    for (int j=1; j<=jmax; j++)
    B(i-1,l)-=D(i-1,j)*B(i+j-1,l);
  }
  }
  return Rcpp::List::create(Rcpp::Named("D")=D,Rcpp::Named("x")=B);
}