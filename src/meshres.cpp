#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP meshres(SEXP vb_, SEXP it_) {
  try {
    NumericMatrix vb(vb_);
    IntegerMatrix it(it_);
    int nit = it.ncol();
    mat vbA(vb.begin(),vb.nrow(),vb.ncol());
    imat itA(it.begin(),it.nrow(),it.ncol());
    vec tmp(3);
    double res = 0.0;
    for (int i=0; i < nit;++i) {
      tmp = vbA.col(it(0,i))-vbA.col(it(1,i));
      res += sqrt(dot(tmp,tmp));
      tmp = vbA.col(it(0,i))-vbA.col(it(2,i));
      res += sqrt(dot(tmp,tmp));
      tmp = vbA.col(it(1,i))-vbA.col(it(2,i));
      res += sqrt(dot(tmp,tmp));
    }
    res /= nit*3;
    return wrap(res);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
