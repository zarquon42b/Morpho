#include "doozers.h"
using namespace Rcpp;

RcppExport SEXP edgePlane(SEXP vb_, SEXP diff_, SEXP edges_) {
  try {
    IntegerMatrix edges(edges_);
    NumericMatrix vb(vb_);
    NumericMatrix diff(diff_);
    unsigned int nedges = edges.nrow();
    mat out(nedges,3); out.zeros();
    std::vector<unsigned int> test;
    for (unsigned int i = 0; i < nedges; i++) {
      int i1 = edges(i,0);
      int i2 = edges(i,1);
      vec vb2 =  vb(i2,_);
      vec diff0 = diff(i2,_);
      double ancath = sqrt(dot(diff0,diff0));
    
      vec resvec = vb(i1,_)-vb(i2,_);
      double angle = angcalcArma(diff0,resvec);
      double lres = sqrt(dot(resvec,resvec));
      resvec = resvec/lres;
      double hypoth = ancath/cos(angle);
      if (hypoth <= lres && hypoth >= 0) {
	out.row(i) = conv_to<rowvec>::from(vb2+hypoth*resvec);
	test.push_back(i);
      }
    }
    uvec myinds = conv_to<uvec>::from(test);
    out = out.rows(myinds);
    return wrap(out);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
  
