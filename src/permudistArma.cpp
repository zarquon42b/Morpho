#include "permudistArma.h"


SEXP permudistArma(SEXP data_, SEXP groups_, SEXP rounds_) {
  try {
    mat armaData = as<mat>(data_);
    arma::ivec armaGroups = Rcpp::as<arma::ivec>(groups_);
    int rounds = Rcpp::as<int>(rounds_);
        
    ivec permuvec = armaGroups;
    int maxlev = armaGroups.max();
    int alldist=0;
    for (int i=1; i < maxlev; ++i)
      alldist +=i;
    List out(alldist);  
    for (int i=0; i < alldist; ++i) {
      NumericVector dist0(rounds+1);
      out[i] =dist0;
    }
    for (int i=0; i <= rounds; ++i) {
      int count = 0;
      if (i > 0)
	permuvec = shuffle(permuvec);
      for (int j0 = 1; j0 < maxlev; ++j0) {
	mat tmp1 = armaData.rows(arma::find(permuvec == j0 ));
	mat mean1 = mean(tmp1,0);
	for(int j1 =j0+1; j1 <= maxlev; ++j1) {
	  mat tmp2 = armaData.rows(arma::find(permuvec == j1 ));
	  mat mean2 = mean(tmp2,0);
	  mat diff = mean1-mean2;
	  double tmpdist = norm(diff,2);
	  NumericVector dists = out[count];
	  dists[i] = tmpdist;
	  out[count]=dists;
	  count +=1;
	}
      }
    }
    return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
