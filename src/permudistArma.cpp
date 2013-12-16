#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP permudistArma(SEXP datar, SEXP groupsr, SEXP roundr, SEXP maxlevr, SEXP alldistr) {
  Rcpp::NumericMatrix data(datar);
  Rcpp::IntegerVector groups(groupsr);
  int rounds = Rcpp::as<int>(roundr);
  int maxlev = Rcpp::as<int>(maxlevr);
  int alldist=0;
  for (int i=1; i < maxlev; ++i)
    alldist +=i;
  int n = data.nrow();
  int m = data.ncol();
  mat armaData(data.begin(), n,m);
  ivec armaGroups(groups.begin(),groups.size(),false);
  ivec permuvec = armaGroups;
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
	double tmpdist = sqrt(dot(diff,diff));
      	NumericVector dists = out[count];
	dists[i] = tmpdist;
        out[count]=dists;
	count +=1;
      }
    }
  }
  return out;
}
