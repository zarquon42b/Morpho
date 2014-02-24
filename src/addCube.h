#ifndef ADDCUBE_H_
#define ADDCUBE_H_
#include <RcppArmadillo.h>
#include "doozers.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

namespace Morpho
{
  template <class mytype>
    class IOCube 
    {
    public:
      static Mat<mytype> addCube(Cube<mytype> mycube) {
	int n = mycube.n_slices;
	Mat<mytype> out = mycube.slice(0);
	if (n > 1) {
	  for (int i=1; i < n; i++) {
	    out += mycube.slice(i);
	  }
	}
	return out;
      }
    };
}
#endif /*ADDCUBE_H_*/
