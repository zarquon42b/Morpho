#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
//using namespace arma;

RcppExport SEXP barycenterCpp(SEXP vb_, SEXP it_) {
  try {
    NumericMatrix vb(vb_);
    IntegerMatrix it(it_);
    int nit = it.ncol();
    NumericMatrix bary(nit,3);
    for (int i=0; i < nit; ++i) {
      bary(i,_) = (vb(_,it(0,i))+vb(_,it(1,i))+vb(_,it(2,i)))/3;
    }
    return wrap(bary);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
