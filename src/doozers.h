#ifndef DOOZERS_H_
#define DOOZERS_H_


#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double angcalcArma(colvec a, colvec b, bool circle = true);

double angcalcRcpp(NumericVector a_, NumericVector b_, bool circle = true);

void crosspArma(colvec x, colvec y, colvec& z);
#endif /*DOOZERS_H_*/
