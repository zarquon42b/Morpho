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
    return 1;
  eig_gen(eigval,X);
  eigval = log(eigval);
  tmp = sqrt(dot(eigval,eigval));
  cdist = real(tmp(0));
  return(cdist);
}
  
RcppExport SEXP covWrap(SEXP s1_, SEXP s2_) {
  NumericMatrix s1(s1_);
  NumericMatrix s2(s2_);
  mat S1(s1.begin(), s1.nrow(), s1.ncol());
  mat S2(s2.begin(), s2.nrow(), s2.ncol());
  double cdist = covDist(S1, S2);
  return wrap(cdist);

}  
