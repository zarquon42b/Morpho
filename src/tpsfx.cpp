#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP tpsfx(SEXP refmat_, SEXP M_, SEXP coefs_, SEXP threads_ = wrap(1)) {
  try {
    typedef unsigned int uint;
    // reference coordinates
    mat refmat = as<mat>(refmat_);
    //M contains homogenous coordinates
    mat M = as<mat>(M_);
    uint lmdim = M.n_cols-1;
    uvec select(lmdim);
    for (uint i = 0; i < lmdim; i++)
      select[i] = i+1;
    // remove leading 1 from homogenous coordinates
    mat Mnohom = M.cols(select);
    mat coefs = as<mat>(coefs_);
    uint m = refmat.n_rows;
    int threads = as<int>(threads_);
    // remove affine coefficients
    mat coefsNoAff = coefs.cols(0, m-1);
    mat result(M.n_rows,coefs.n_rows); result.zeros();
    
    #pragma omp parallel for schedule(static) num_threads(threads)

    for (uint i=0; i < Mnohom.n_rows; ++i) {
      colvec x(m);
      for (uint j=0; j < m; ++j) {
	mat tmp = refmat.row(j) - Mnohom.row(i);
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
      vec tmpres = coefs.cols(m,m+lmdim)*M.row(i).t();
      result.row(i) = (tmp+tmpres).t();
    }
    return wrap(result);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
