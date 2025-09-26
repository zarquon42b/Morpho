#ifndef angcal_H_
#define angcal_H_

#include <RcppArmadillo.h>
#include "doozers.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

double angcalcArma(vec a, vec b) {
  try {
    double alen = norm(a,2);
    double blen = norm(b,2);
    if (alen > 0)
      a = a/alen;
    if (blen > 0)
      b = b/blen;
    vec diffvec = a-b;
    double angle = acos((dot(diffvec,diffvec)-2)/-2);
    return angle;
  } catch (std::exception& e) {
    forward_exception_to_r( e );
  } catch (...) {
    ::Rf_error("unknown exception");
  } return 1.0; // -Wall
}
		
double angcalcRcpp(NumericVector a_, NumericVector b_) {
  try {
    colvec a(a_.begin(),a_.size(),false);
    colvec b(b_.begin(),b_.size(),false);
    double alen = sqrt(dot(a,a));
    double blen = sqrt(dot(b,b));
    if (alen > 0)
      a = a/alen;
    if (blen > 0)
      b = b/blen;
    colvec diffvec = a-b;
    double angle = acos((dot(diffvec,diffvec)-2)/-2);
    return angle;
  } catch (std::exception& e) {
    forward_exception_to_r( e );
  } catch (...) {
    ::Rf_error("unknown exception");
  } return 1.0; 
}
void crosspArma(colvec x, colvec y, colvec& z) {
  z(0) = x(1)*y(2)-x(2)*y(1);
  z(1) = x(2)*y(0)-x(0)*y(2);
  z(2) = x(0)*y(1)-x(1)*y(0);
  double lz = sqrt(dot(z,z));
  if (lz > 0) 
    z = z/lz;
}


double tpskernel(mat tmp, int dim) {
  double out;
  if(dim > 2)
    double out = -sqrt(dot(tmp,tmp));
  else {
	  double tmp0 = dot(tmp,tmp);
	  if (tmp0 == 0) 
	    out = 0;
	  else 
	    out = tmp0*log(tmp0);
	}
  return(out);
}


//2d h^4*log(h);
//3d h^3
double tpskernelCube(mat tmp, int dim) {
  double out;
  double h = sqrt(dot(tmp,tmp));
  if(dim > 2)
    double out = pow(h,3.0);
  else {
    double tmp0 = pow(h,4.0);
	  if (tmp0 == 0) 
	    out = 0;
	  else 
	    out = tmp0*log(h);
	}
  return(out);
}


#endif /*angcal_H_*/
