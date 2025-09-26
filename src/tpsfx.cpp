#include <RcppArmadillo.h>
#include "doozers.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace std;
using namespace arma;




RcppExport SEXP tpsfx(SEXP refmat_, SEXP M_, SEXP coefs_, SEXP tpskernel_, SEXP threads_ = wrap(1)) {
  try {
    typedef unsigned int uint;
    // reference coordinates
    mat refmat = as<mat>(refmat_);
    //M contains homogenous coordinates
    mat M = as<mat>(M_);
    uint lmdim = M.n_cols-1;
    int tpskerneltype = as<int>(tpskernel_);
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
    Rprintf("%i\n",tpskerneltype);
    #pragma omp parallel for schedule(static) num_threads(threads)

    for (uint i=0; i < Mnohom.n_rows; ++i) {
      colvec x(m);
      for (uint j=0; j < m; ++j) {
	mat tmp = refmat.row(j) - Mnohom.row(i);
	if (tpskerneltype == 0)
	  x(j) = tpskernel(tmp,lmdim);
	else {
	  x(j) = tpskernelCube(tmp,lmdim);
	}
      }
      vec tmp = coefsNoAff*x;
      vec tmpres = coefs.cols(m,m+lmdim)*M.row(i).t();
      result.row(i) = (tmp+tmpres).t();
    }
    return wrap(result);
  } catch (std::exception& e) {
    forward_exception_to_r( e );
  } catch (...) {
    ::Rf_error("unknown exception");
  } return R_NilValue; 
}

