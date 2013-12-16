#include "angcalc.h"
using namespace Rcpp;

RcppExport SEXP ang_calcC(SEXP x_, SEXP y_, SEXP circle_) {
  NumericVector x(x_);
  NumericVector y(y_);
  bool circle = as<bool>(circle_);
  colvec xA(x.begin(), x.size());
  colvec yA(y.begin(), y.size());
  double pi = 3.141592653589793239;
  double rho = angcalcRcpp(x, y);
  if (! circle)
    rho = pi -rho;
  return wrap(rho);
}

   
  
