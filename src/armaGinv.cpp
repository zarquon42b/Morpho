#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;


RcppExport SEXP armaGinv(SEXP matIn_, double *tol) {
  NumericMatrix matIn(matIn_);
  mat matA(matIn.begin(), matIn.nrow(), matIn.ncol());
  mat invA = pinv(matA, *tol);
  return wrap(invA);
  

}
