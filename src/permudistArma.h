#ifndef PERMUDIST_H_
#define PERMUDIST_H_

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP permudistArma(SEXP data_, SEXP groups_, SEXP rounds_);

#endif /*PERMUDIST_H_*/
