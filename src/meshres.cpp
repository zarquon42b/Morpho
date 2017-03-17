#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP meshresCpp(SEXP vb_, SEXP it_) {
  try {
        
    arma::mat vbA = Rcpp::as<arma::mat>(vb_);
    //mat vbA(vb.begin(),vb.nrow(),vb.ncol());
    imat itA = as<arma::imat>(it_);
    int nit = itA.n_cols;
    //imat itA(it.begin(),it.nrow(),it.ncol());
    vec tmp(3);
    double res = 0.0;
    for (int i=0; i < nit;++i) {
      tmp = vbA.col(itA(0,i))-vbA.col(itA(1,i));
      res += sqrt(dot(tmp,tmp));
      tmp = vbA.col(itA(0,i))-vbA.col(itA(2,i));
      res += sqrt(dot(tmp,tmp));
      tmp = vbA.col(itA(1,i))-vbA.col(itA(2,i));
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
