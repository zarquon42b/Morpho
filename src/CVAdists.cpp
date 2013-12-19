#include "CVAdists.h"

SEXP CVAdists(SEXP data_, SEXP  groups_, SEXP rounds_, SEXP winv_) {
  NumericMatrix data(data_);
  NumericMatrix winv(winv_);
  IntegerVector groups(groups_);
  int rounds = as<int>(rounds_);
  //bool n3 = as<bool>(n3_);
  int n = data.nrow();
  int m = data.ncol();
  mat armaData(data.begin(), n,m);
  mat winvA(winv.begin(), winv.nrow(), winv.ncol());
  ivec armaGroups(groups.begin(),groups.size(),false);
  int maxlev = armaGroups.max();
  int alldist=0;
  for (int i=1; i < maxlev; ++i)
    alldist +=i;
  
  ivec permuvec = armaGroups;
  List outPlain(alldist);
  List outMaha(alldist);
  //setup output lists and fill with empty vectors
  for (int i=0; i < alldist; ++i) {
    NumericVector dist0(rounds+1);
   outPlain[i] =dist0;
   NumericVector dist1(rounds+1);
   outMaha[i] = dist1;
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
      	NumericVector dists = outPlain[count];
	dists[i] = tmpdist;
	outPlain[count] = dists;
	// mahalanobis distances
	mat tmpdist0 = sqrt(diff*winvA*diff.t());
	dists = outMaha[count];
	dists[i] = tmpdist0(0,0);
        outMaha[count] = dists;
	count +=1;
      }
    }
  }
 return List::create(Named("Maha")=outMaha,
		     Named("Plain") = outPlain)
   ;
}
