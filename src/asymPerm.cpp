#include "asymPerm.h"
#include "doozers.h"


using namespace Rcpp;
using namespace std;
using namespace arma;

SEXP asymPerm(SEXP asymr, SEXP groupsr, SEXP roundr) {
  try {
    mat armaAsym =as<mat>(asymr);
    arma::ivec armaGroups = Rcpp::as<arma::ivec>(groupsr);
    int rounds = Rcpp::as<int>(roundr);
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
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP asymPermuteCpp(SEXP data_, SEXP groups_, SEXP rounds_) {
  try {
    mat armaData = as<mat>(data_);
    arma::ivec armaGroups = Rcpp::as<arma::ivec>(groups_);
    int rounds = Rcpp::as<int>(rounds_);
    ivec permuvec = armaGroups;
    int maxlev = armaGroups.max();
    int alldist=0;
    for (int i=1; i < maxlev; ++i)
      alldist +=i;
    List outdiff(alldist);  
    List angdiff(alldist);
    for (int i=0; i < alldist; ++i) {
      NumericVector dist0(rounds+1);
      NumericVector dist1(rounds+1);
      outdiff[i] = dist0;
      angdiff[i] = dist1;
    }
    for (int i=0; i <= rounds; ++i) {
      int count = 0;
      if (i > 0)
	permuvec = shuffle(permuvec);
      for (int j0 = 1; j0 < maxlev; ++j0) {
	mat tmp1 = armaData.rows(arma::find(permuvec == j0 ));
	mat mean1mat = mean(tmp1,0);
	vec mean1 = vectorise(mean1mat);
	double mean1len = sqrt(dot(mean1,mean1));
	for(int j1 =j0+1; j1 <= maxlev; ++j1) {
	  mat tmp2 = armaData.rows(arma::find(permuvec == j1 ));
	  mat mean2mat = mean(tmp2,0);
	  vec mean2 = vectorise(mean2mat);
	  double mean2len = sqrt(dot(mean2,mean2));
	  double ang = angcalcArma(mean1,mean2);
	  double tmpdist = abs(mean1len-mean2len);
	  NumericVector dists = outdiff[count];
	  NumericVector angs = angdiff[count];
	  dists[i] = tmpdist;
	  angs[i] = ang;
	  outdiff[count]=dists;
	  angdiff[count] = angs;
	  count +=1;
	}
      }
    }
    return List::create(Named("angles") = angdiff,
			Named("dists") = outdiff
			);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

