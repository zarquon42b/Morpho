#ifndef updateNormals_H_
#define updateNormals_H_

#include "doozers.h"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP  updateVertexNormals(SEXP vb_, SEXP it_,SEXP angweight_);

RcppExport SEXP updateFaceNormals(SEXP vb_, SEXP it_);


#endif /*updateNormals_H_*/
