#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP createL(SEXP Matrix_, SEXP threads_= wrap(1)) {
  try {
    mat MatrixA = as<mat>(Matrix_);
    int threads = as<int>(threads_);
    int k = MatrixA.n_rows;
    mat K(k,k); K.zeros();
    int m = MatrixA.n_cols;
    int omp = 1;
    if (threads == 1)
      omp = 0;
#pragma omp parallel for schedule(static) num_threads(threads) if (omp)
    for (int i=0; i < (k-1); ++i) {
      for(int j=(i+1); j < k; ++j) {
	mat diff = MatrixA.row(i)-MatrixA.row(j);
	if (m == 3)
	  K(i,j) = -sqrt(dot(diff,diff));
	if (m == 2) {
	  double r2 = dot(diff,diff);
	  K(i,j) = r2*log(r2);
	}
      }
    }
    K = K+K.t();
    return wrap(K);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
  
