#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

double check_inner(std::vector<double> poly_x, std::vector<double> poly_y, double a, double b, double x0, double y0);

vec get_jumperpoint(std::vector<double> poly_x, std::vector<double> poly_y, double sizeup,double a, double b, double x0,double y0);

RcppExport SEXP inscribeEllipseCpp(SEXP polyMat_, SEXP step_, SEXP iters_,SEXP init_point_);
