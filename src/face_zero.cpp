#include <Rcpp.h> 

using namespace Rcpp;
using namespace std;

RcppExport SEXP face_zero(SEXP it_) {
  try {
    IntegerMatrix it(it_);
    int nit = it.ncol();
    IntegerVector out(nit);
    for (int i=0; i < nit; i++) {
       out(i) = 1;
      for(int j=0; j < it.nrow(); j++) {
	if (it(j,i) == 0)
	  out(i) = 0;
	
      }
    }
      //IntegerVector tmp = it(_,i);
      //out(i) = std::accumulate(tmp.begin(), tmp.end(), 1,std::multiplies<int>());
      //Rprintf("%i\n",a);
      //out(i) = a;
    
    return wrap(out);
  } catch (std::exception& e) {
    forward_exception_to_r( e );
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
