#ifndef CVA_DISTS_H_
#define CVA_DISTS_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

RcppExport SEXP CVAdists(SEXP data_, SEXP  groups_, SEXP rounds_, SEXP winv_);

#endif /*CVA_DISTS_H_*/
