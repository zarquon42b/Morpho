#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

RcppExport SEXP createL(SEXP Matrix_) {
   NumericMatrix Matrix(Matrix_);
   int m = Matrix.ncol();
   int k = Matrix.nrow();
   mat MatrixA(Matrix.begin(), Matrix.nrow(), Matrix.ncol());
   mat K(k,k);
   for (int i=0; i < k; ++i) {
     for(int j=0; j < k; ++j) {
       mat diff = MatrixA.row(i)-MatrixA.row(j);
       K(i,j) = -sqrt(dot(diff,diff));
     }
   }
   return wrap(K);
}
  
