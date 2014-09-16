#include <Rcpp.h> 

using namespace Rcpp;
using namespace std;

RcppExport SEXP face_zero(SEXP it_) {
  try {
    IntegerMatrix it(it_);
    int nit = it.ncol();
    IntegerVector out(nit);
    for (int i=0; i < nit; i++) {
      if (it(0,i) == 0 || it(1,i) == 0 || it(2,i) == 0)
	out(i) = 0;
      else 
	out(i) = 1;
      //IntegerVector tmp = it(_,i);
      //out(i) = std::accumulate(tmp.begin(), tmp.end(), 1,std::multiplies<int>());
      //Rprintf("%i\n",a);
      //out(i) = a;
    }
    return wrap(out);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
