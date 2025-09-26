#include <RcppArmadillo.h>
#include "doozers.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP createL(SEXP Matrix_,SEXP tpskernel_, SEXP threads_= wrap(1)) {
  try {
    mat MatrixA = as<mat>(Matrix_);
    int threads = as<int>(threads_);
    int tpskerneltype = as<int>(tpskernel_);
    int k = MatrixA.n_rows;
    mat K(k,k); K.zeros();
    int m = MatrixA.n_cols;
    int omp = 1;
    Rprintf("CreateL %i\n",tpskerneltype);
    Rprintf("CreateL %i\n",m);
    if (threads == 1)
      omp = 0;
#pragma omp parallel for schedule(static) num_threads(threads) if (omp)
    for (int i=0; i < (k-1); ++i) {
      for(int j=(i+1); j < k; ++j) {
	mat diff = MatrixA.row(i)-MatrixA.row(j);
	if (tpskerneltype == 0)
	  K(i,j) = tpskernel(diff,m);
	else {
	  K(i,j) = tpskernelCube(diff,m);
	}
      }
    }
    K = K+K.t();
    return wrap(K);
  } catch (std::exception& e) {
    forward_exception_to_r( e );
  } catch (...) {
    ::Rf_error("unknown exception");
  } return R_NilValue; 
}
  
