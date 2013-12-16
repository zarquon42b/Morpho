#ifndef updateNormals_H_
#define updateNormals_H_

#include "angcalc.h"

RcppExport SEXP  updateNormals(SEXP vb_, SEXP it_,SEXP angweight_);

RcppExport SEXP updateFaceNormals(SEXP vb_, SEXP it_);


#endif /*updateNormals_H_*/
