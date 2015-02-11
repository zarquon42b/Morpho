#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP tpsfx(SEXP A_,SEXP B_,SEXP Bh_, SEXP coefs_) {
  try {
    typedef unsigned int uint;
    NumericMatrix A(A_);
    NumericMatrix B(B_);
    NumericMatrix Bh(Bh_);
    NumericMatrix coefs(coefs_);
    uint m = A.nrow();
    uint lmdim = A.ncol();
    mat AA(A.begin(), A.nrow(), A.ncol());
    mat BA(B.begin(), B.nrow(), B.ncol());
    mat BhA(Bh.begin(), Bh.nrow(), Bh.ncol());
    mat coefsA(coefs.begin(), coefs.nrow(), coefs.ncol());
    mat coefsNoAff = coefsA.cols(0, m-1);
    mat result = BA; result.zeros();
    colvec x(m);
  
    for (uint i=0; i < BA.n_rows; ++i) { 
      for (uint j=0; j < m; ++j) {
	mat tmp = AA.row(j) - BA.row(i);
	if (lmdim > 2) {
	  x(j) = -sqrt(dot(tmp,tmp));
	} else {
	  double tmp0 = dot(tmp,tmp);
	  if (tmp0 == 0) 
	    x(j) = 0;
	  else 
	    x(j) = tmp0*log(tmp0);
	}
      }
      vec tmp = coefsNoAff*x;
      vec tmpres = coefsA.cols(m,m+lmdim)*BhA.row(i).t();
      result.row(i) = (tmp+tmpres).t();
    }
    return wrap(result);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
