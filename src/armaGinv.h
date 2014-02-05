#ifndef _armaGinv_H
#define _armaGinv_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP armaGinv(SEXP matIn_, SEXP tol_);

#endif
