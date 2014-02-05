#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;


double covDist(mat s1, mat s2) {
  double cdist;
  mat X, EigVec;
  cx_vec eigval, tmp(1);
  bool check = solve(X, s1, s2);
  if (!check)
    return 2.0e30;
  eig_gen(eigval,X);
  eigval = log(eigval);
  tmp = sqrt(dot(eigval,eigval));
  cdist = real(tmp(0));
  return(cdist);
}

mat covPCA(mat data, ivec groups, int rounds, bool scramble) {
  typedef unsigned int uint;
  uint maxlev = groups.max();
  mat dists(maxlev,maxlev);
  double check;
  List covaList(maxlev);
 
    for (uint i = 0; i < maxlev; ++i) {
      if (!scramble) {
      covaList[i] = cov(data.rows(arma::find(groups == (i+1))));
      } else {

	mat tmpdat = data.rows(arma::find(groups == (i+1)));
	uint nrow = tmpdat.n_rows;
	uvec shaker = randi<uvec>(nrow, distr_param(0,nrow));
	covaList[i] = cov(tmpdat.rows(shaker));
      }
    }
  dists.zeros();
  for (uint i = 0; i < (maxlev-1); ++i) {
    for (uint j = i+1; j < (maxlev); ++j){
      mat tmp0 = covaList[i];
      mat tmp1 = covaList[j];
      check = covDist(tmp0,tmp1);
      dists(j,i) = check;
    }
   
  }
 vec checkerr = conv_to<vec>::from(dists);
 if (any(checkerr == 2.0e30)) {
   mat errout(0,0);
   return errout;
 } else
   return dists;
  
}
cube covPCAboot(mat data, ivec groups, int rounds) {
  
  typedef unsigned int uint;
  uint maxlev = groups.max();  
  cube alldist(maxlev, maxlev, rounds);
  for (int i = 0; i < rounds;){
    mat result = covPCA(data, groups, 0, true);
    if (result.n_cols > 0) {
      alldist.slice(i) = result;
      i++;//only increment if covPCA did not fail
    }
  }

    return alldist;

}


    
RcppExport SEXP covPCAwrap(SEXP data_, SEXP groups_, SEXP scramble_) {
  //bool scramble = Rcpp::as<bool>(scramble_);
 int scramble = Rcpp::as<int>(scramble_);
  Rcpp::NumericMatrix data(data_);
  Rcpp::IntegerVector groups(groups_);
  mat armaData(data.begin(), data.nrow(),data.ncol());
  ivec armaGroups(groups.begin(),groups.size(),false);
  cube out = covPCAboot(armaData,armaGroups,scramble);
  return wrap(out);

}
  
RcppExport SEXP covWrap(SEXP s1_, SEXP s2_) {
  NumericMatrix s1(s1_);
  NumericMatrix s2(s2_);
  mat S1(s1.begin(), s1.nrow(), s1.ncol());
  mat S2(s2.begin(), s2.nrow(), s2.ncol());
  double cdist = covDist(S1, S2);
  return wrap(cdist);

}  
