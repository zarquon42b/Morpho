#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

RcppExport SEXP fastSubsetMeans(SEXP x_, SEXP inds_,SEXP threads_) {
   try {
     mat x = as<mat>(x_);
  uvec inds = as<uvec>(inds_);
  int threads = as<int>(threads_);
  int maxlev = inds.max();
  mat center(maxlev,x.n_cols);
#pragma omp parallel for schedule(static) num_threads(threads)
  for (int i =0; i < maxlev;i++) {
    
    uvec tmpinds = find(inds ==i);
    mat tmpmat = x.rows(tmpinds);
    rowvec tmpresult(x.n_cols);tmpresult.fill(0);
    for (int j = 0; j < tmpinds.size();j++)
      tmpresult += tmpmat.row(j);
    tmpresult /= tmpinds.size();
    center.row(i) = tmpresult;
  }
  return wrap(center);
   }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
