#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP meanMat(SEXP matIn_) {
  try {  
    NumericMatrix matIn(matIn_);
    mat matA(matIn.begin(), matIn.nrow(), matIn.ncol());
    mat mout = mean(matA);
    return wrap(mout);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
