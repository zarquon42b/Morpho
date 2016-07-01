#include "doozers.h"
using namespace Rcpp;

RcppExport SEXP ang_calcC(SEXP x_, SEXP y_) {
  try{
    NumericMatrix x(x_);
    NumericVector angle(x.nrow());
    NumericVector y(y_);
    for (int i = 0; i < x.nrow(); i++) {
      angle(i) = angcalcRcpp(x(i,_),y);
    }

    return wrap(angle);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

   
RcppExport SEXP ang_calcM(SEXP x_, SEXP y_) {
  try{
    NumericMatrix x(x_);
    NumericVector angle(x.nrow());
    NumericMatrix y(y_);
    
    for (int i = 0; i < x.nrow(); i++) {
      angle(i) = angcalcRcpp(x(i,_),y(i,_));
    }

    return wrap(angle);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

