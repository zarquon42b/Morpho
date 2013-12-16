#include "asymPerm.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

SEXP asymPerm(SEXP asymr, SEXP groupsr, SEXP roundr) {
Rcpp::NumericMatrix asym(asymr);
  Rcpp::IntegerVector groups(groupsr);
  int rounds = Rcpp::as<int>(roundr);
  int n = asym.nrow();
  int m = asym.ncol();
  mat armaAsym(asym.begin(), n,m);
  ivec armaGroups(groups.begin(),groups.size(),false);
  ivec permuvec = armaGroups;
  std::vector<double> diff,angle;
  for (int i=0; i <= rounds; ++i) {
    if (i > 0)
        permuvec = shuffle(permuvec);
    mat tmp1 = armaAsym.rows(arma::find(permuvec == 1 ));
    mat mean1 = mean(tmp1,0);
    double mean1len = sqrt(dot(mean1,mean1));
    mat tmp2 = armaAsym.rows(arma::find(permuvec == 2 ));
    mat mean2 = mean(tmp2,0);
    double mean2len = sqrt(dot(mean2,mean2));
    mean1 = mean1/mean1len;
    mean2 = mean2/mean2len;
    mat meandiff = mean1-mean2;
    double ang = acos((dot(meandiff,meandiff)-2)/-2);
    angle.push_back(ang);
    diff.push_back(abs(mean1len-mean2len));
      }
  return List::create(Named("diff")=diff,
		      Named("angle")=angle
		      );
}
