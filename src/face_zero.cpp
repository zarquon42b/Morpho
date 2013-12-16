#include <Rcpp.h> 

using namespace Rcpp;
using namespace std;

RcppExport SEXP face_zero(SEXP it_) {
  IntegerMatrix it(it_);
  int nit = it.ncol();
  IntegerVector out(nit);
  for (int i=0; i < nit; ++i) {
    out(i) = it(0,i)*it(1,i)*it(2,i);
  }
  return wrap(out);
}
