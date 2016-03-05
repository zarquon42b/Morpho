#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP tpsfx(SEXP A_,SEXP B_,SEXP Bh_, SEXP coefs_, SEXP threads_ = wrap(1)) {
  try {
    typedef unsigned int uint;
    mat AA = as<mat>(A_);
    mat BA = as<mat>(B_);
    mat BhA = as<mat>(Bh_);
    mat coefsA = as<mat>(coefs_);
    NumericMatrix A(A_);
    NumericMatrix B(B_);
    NumericMatrix Bh(Bh_);
    NumericMatrix coefs(coefs_);
    uint m = AA.n_rows;
    uint lmdim = AA.n_cols;
    int threads = as<int>(threads_);
    
    mat coefsNoAff = coefsA.cols(0, m-1);
    mat result = BA; result.zeros();
    
    uint i;
    #pragma omp parallel for schedule(static) private(i) num_threads(threads)

    for (i=0; i < BA.n_rows; ++i) {
      colvec x(m);
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
