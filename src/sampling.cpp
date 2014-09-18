#include "sampling.h"

using namespace Rcpp;
using namespace arma;

arma::uvec  myrandu(unsigned int n, unsigned int from, unsigned int to) {
  NumericVector out = floor(runif(n,from,to));
  vec out1(out.begin(),out.size(),false);
  uvec out2 = conv_to<uvec>::from(out1);
  return out2;
}
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

arma::ivec randomShuffle(arma::ivec b) {
    // already added by sourceCpp(), but needed standalone
    Rcpp::RNGScope scope;             

    std::random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}
