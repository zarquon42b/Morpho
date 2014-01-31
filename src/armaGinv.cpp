#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;


RcppExport SEXP armaGinv(SEXP matIn_, SEXP tol_) {
  if (!Rf_isMatrix(matIn_)){
    return wrap(1);
  } else {
  NumericMatrix matIn(matIn_);
  mat matA(matIn.begin(), matIn.nrow(), matIn.ncol());
  mat invA;
  bool check;
  if (Rf_isNumeric(tol_)) {
    double tol = as<double>(tol_);
    check =pinv(invA, matA, tol);
  } else {
    check = pinv(invA, matA);
  }
  if (check)
    return wrap(invA);
  else 
    return wrap(1);
  }

}
