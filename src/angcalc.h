#ifndef angcal_H_
#define angcal_H_


#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

double angcalcArma(colvec a, colvec b);
		
double angcalcRcpp(NumericVector a_, NumericVector b_);

void crosspArma(colvec x, colvec y, colvec& z);
#endif /*angcal_H_*/
