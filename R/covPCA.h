#ifndef _covPCA_H
#define _covPCA_H
#include <RcppArmadillo.h>


double covDist(mat &s1, mat &s2);

mat covDistMulti(mat data, ivec groups, bool scramble);

cube covPCAboot(mat data, ivec groups, int rounds);

cube covPCApermute(mat data, ivec groups, int rounds);

List covMDS(mat &dists);

RcppExport SEXP covPCAwrap(SEXP data_, SEXP groups_, SEXP scramble_, SEXP rounds_);

RcppExport SEXP covWrap(SEXP s1_, SEXP s2_) {
  NumericMatrix s1(s1_);
  NumericMatrix s2(s2_);
  mat S1(s1.begin(), s1.nrow(), s1.ncol());
  mat S2(s2.begin(), s2.nrow(), s2.ncol());
  double cdist = covDist(S1, S2);
  return wrap(cdist);

}  
#endif
