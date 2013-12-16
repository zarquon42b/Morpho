#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP tpsfx(SEXP A_,SEXP B_,SEXP Bh_, SEXP coefs_) {
  NumericMatrix A(A_);
  NumericMatrix B(B_);
  NumericMatrix Bh(Bh_);
  NumericMatrix coefs(coefs_);
  int m = A.nrow();
  mat AA(A.begin(), A.nrow(), A.ncol());
  mat BA(B.begin(), B.nrow(), B.ncol());
  mat BhA(Bh.begin(), Bh.nrow(), Bh.ncol());
  mat coefsA(coefs.begin(), coefs.nrow(), coefs.ncol());
  mat coefsNoAff = coefsA.cols(0, m-1);
  mat result = BA; result.zeros();
  colvec x(m);
  
  for (int i=0; i < BA.n_rows; ++i) { 
    for (int j=0; j < m; ++j) {
      mat tmp = AA.row(j) - BA.row(i);
      x(j) = -sqrt(dot(tmp,tmp));
    }
    vec tmp = coefsNoAff*x;
    vec tmpres = coefsA.cols(m,m+3)*BhA.row(i).t();
    result.row(i) = (tmp+tmpres).t();
  }
  return wrap(result);
}
