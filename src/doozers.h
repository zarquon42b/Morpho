#ifndef DOOZERS_H_
#define DOOZERS_H_
#ifndef ARMA_DONT_PRINT_ERRORS 
#define ARMA_DONT_PRINT_ERRORS 
#endif

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double angcalcArma(colvec a, colvec b);
		
double angcalcRcpp(NumericVector a_, NumericVector b_);

void crosspArma(colvec x, colvec y, colvec& z);

double tpskernel(mat tmp, int dim);

double tpskernelCube(mat tmp, int dim);

#endif /*DOOZERS_H_*/
