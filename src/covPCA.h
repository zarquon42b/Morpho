#ifndef _covPCA_H
#define _covPCA_H
#ifndef ARMA_DONT_PRINT_ERRORS 
#define ARMA_DONT_PRINT_ERRORS 
#include <RcppArmadillo.h>
#endif

using namespace Rcpp;
using namespace std;
using namespace arma;

double covDist(mat &s1, mat &s2);

mat covDistMulti(mat &data, ivec groups, bool scramble);

cube covPCAboot(mat &data, ivec groups, int rounds);

cube covPCApermute(mat &data, ivec groups, int rounds);

List covMDS(mat &dists);

RcppExport SEXP covPCAwrap(SEXP data_, SEXP groups_, SEXP scramble_, SEXP rounds_);

RcppExport SEXP covWrap(SEXP s1_, SEXP s2_);   

#endif
