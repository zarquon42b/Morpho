#include "armaGinvCpp.h"

SEXP armaGinvCpp(SEXP matIn_, SEXP tol_) {
  try {
    if (!Rf_isMatrix(matIn_)){
      return wrap(1);
    } else {
      mat matA = as<mat>(matIn_);
      mat invA;
      bool check;
      if (Rf_isNumeric(tol_)) {
	double tol = as<double>(tol_);
	check =pinv(invA, matA, tol);
      } else {
	check = pinv(invA, matA);
      }
      if (check)
	return wrap(invA);
      else 
	return wrap(1);
  
    }
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
