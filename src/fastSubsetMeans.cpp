#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

RcppExport SEXP fastSubsetMeans(SEXP x_, SEXP inds_, SEXP k_, SEXP threads_) {
  try {
    mat x = as<mat>(x_);
    int k = as<int>(k_);
    uvec inds = as<uvec>(inds_);
    int threads = as<int>(threads_);
    mat center(k,x.n_cols);
    vec checkempty(k);checkempty.fill(0);
  
    center.fill(0);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int i =0; i < k;i++) {
      uvec tmpinds = find(inds ==i);
      mat tmpmat = x.rows(tmpinds);
      rowvec tmpresult(x.n_cols);tmpresult.fill(0);
      if (tmpinds.size() == 0)
	checkempty(i) = 1;
      for (int j = 0; j < tmpinds.size();j++)
	tmpresult += tmpmat.row(j);
      tmpresult /= tmpinds.size();
      center.row(i) = tmpresult;
    }
    List out = List::create(Named("centers")=center,
			    Named("checkempty")=checkempty);
				
    return out;
  }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
