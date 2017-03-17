#ifndef _armaGinvCpp_H
#define _armaGinvCpp_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP armaGinvCpp(SEXP matIn_, SEXP tol_);

#endif
