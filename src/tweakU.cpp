#include <Rcpp.h>
using namespace Rcpp;
using std::vector;

void testandpush(vector<int>& rows, vector<int>& cols, vector<double>& x, double testit, int r, int c) {
  if (testit != 0) {
    rows.push_back(r);
    cols.push_back(c);
    x.push_back(testit);
  }
}

RcppExport SEXP tweakU(SEXP tanvec_, SEXP m_, SEXP type_, SEXP SMsort_) {
  try {
    NumericMatrix tanvec(tanvec_);
    IntegerVector SMsort(SMsort_);
    int m = as<int>(m_);
    int type = as<int>(type_);
    int k = tanvec.rows();
    std::vector<int> rows, cols;
    std::vector<double> x;
    for (int i = 0; i < m; i++) {
      testandpush(rows,cols,x,tanvec(SMsort[i]-1,0),SMsort[i],i);
      testandpush(rows,cols,x,tanvec(SMsort[i]-1,1),k+SMsort[i],i);
      testandpush(rows,cols,x,tanvec(SMsort[i]-1,2),2*k+SMsort[i],i);
      if (type == 1 || type == 2) {
	testandpush(rows,cols,x,tanvec(SMsort[i]-1,3),SMsort[i],i+m);
	testandpush(rows,cols,x,tanvec(SMsort[i]-1,4),k+SMsort[i],i+m);
	testandpush(rows,cols,x,tanvec(SMsort[i]-1,5),2*k+SMsort[i],i+m);
      }
      if (type == 2) {
	testandpush(rows,cols,x,tanvec(SMsort[i]-1,6),SMsort[i],i+2*m);
	testandpush(rows,cols,x,tanvec(SMsort[i]-1,7),k+SMsort[i],i+2*m);
	testandpush(rows,cols,x,tanvec(SMsort[i]-1,8),2*k+SMsort[i],i+2*m);
      }
    }
    return List::create(Named("rows") = rows,
			Named("cols") = cols,
			Named("x") = x
			);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  } 
}
			
	
      
