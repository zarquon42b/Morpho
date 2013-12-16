#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
//using namespace arma;

RcppExport SEXP barycenter(SEXP vb_, SEXP it_) {
  NumericMatrix vb(vb_);
  IntegerMatrix it(it_);
  int nit = it.ncol();
  NumericMatrix bary(nit,3);
  for (int i=0; i < nit; ++i) {
    bary(i,_) = (vb(_,it(0,i))+vb(_,it(1,i))+vb(_,it(2,i)))/3;
  }
  return wrap(bary);
}
