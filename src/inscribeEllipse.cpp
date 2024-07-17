#include <RcppArmadillo.h>
#include "inscribeEllipse.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

double check_inner(std::vector<double> poly_x, std::vector<double> poly_y, double a, double b, double x0, double y0) {
  
  double minrr = 5;
  for (int i=0; i < poly_x.size();i++) {
    double x = poly_x[i];
    double y = poly_y[i];
    double rr = pow((x-x0)/a,2) + pow((y-y0)/b,2);
    if (rr < minrr)
      minrr = rr;
  }
  return minrr;
}

vec get_jumperpoint(std::vector<double> poly_x, std::vector<double> poly_y, double sizeup,double a, double b, double x0,double y0) {
  std::vector<double> badx ,bady;
    
  for (int i=0; i < poly_x.size();i++) {
    double x = poly_x[i];
    double y = poly_y[i];
    double rr = pow((x-x0)/(a+sizeup),2)+pow(((y-y0)/(b+sizeup)),2);
    if (rr < 1) {
      badx.push_back(x);
      bady.push_back(y);
    }
  }
  vec badxA(badx);
  vec badyA(bady);
  double xmean = mean(badxA);
  double ymean = mean(badyA);
  vec out(2);
  out[0] = x0-xmean;
  out[1] = y0-ymean;
 
  
  return out;
};




RcppExport SEXP inscribeEllipseCpp(SEXP poly_x_, SEXP poly_y_, SEXP step_, SEXP iters_,SEXP init_point_) {
  
  try{
    std::vector<double> px_old = as<std::vector<double>>(poly_x_);
    std::vector<double> py_old = as<std::vector<double>>(poly_y_);
    double step = as<double>(step_);
    int iters = as<int>(iters_);
    std::vector<double> init_point = as<std::vector<double>>(init_point_);
    std::vector<double> bestxy = init_point;
    double init_radius = step;
  
    std::vector<double> px, py;
    for (int i = 0; i < (px_old.size()-1); i++) {
      double dx = px_old[i+1]-px_old[i];
      double dy = py_old[i+1]-py_old[i];
      double len1 = pow((dx*dx+dy*dy),0.5);
      px.push_back(px_old[i]);
      py.push_back(py_old[i]);
    
      if (len1 >= step){
	int count = len1/step;
	for (int ii=0; ii < count; ii++ ) {
	  px.push_back(px_old[i]+ 1.0*ii/count*dx);
	  py.push_back(py_old[i]+ 1.0*ii/count*dy);
	}
      }
    }
    std::vector<double> pxOldRed = px_old;
    pxOldRed.erase(pxOldRed.begin());
    px.insert(px.end(),pxOldRed.begin(),pxOldRed.end());
    std::vector<double> pyOldRed = py_old;
    pyOldRed.erase(pyOldRed.begin());
    py.insert(py.end(),pyOldRed.begin(),pyOldRed.end());
    double minrr = 6;
    double rx = init_radius;
    double besta = init_radius;
    double bestb = init_radius;
    double ry = init_radius;
    double bestx = init_point[0];
    double xc = init_point[0];
    double besty = init_point[1];
    double yc = init_point[1];
    double maxarea = step*step;
    double bestiter = 0;
   
 
    for (int iterat = 0; iterat < iters; iterat++) {
   
      double s1 = step;
      double s3 = step;
      double s2 = step;
      double rxprev = rx;
      double ryprev = ry;
      double xcprev = xc;
      double ycprev = yc;
     

      //double ci =  check_inner(px,py,rx+s1,ry+s2,xc,yc) ;
      //Rprintf("%f\n",ci);
      
      if (check_inner(px,py,rx+s1,ry+s2,xc,yc) >= 1) { 
	rx = rx+s1;
	ry = ry+s2;
      }
     
	if (check_inner(px,py,rx+s1,ry,xc,yc) >= 1) 
	rx = rx+s1;

	if (check_inner(px,py,rx,ry+s1,xc,yc) >= 1)
	ry = ry+s1;

	if (check_inner(px,py,rx+s1,ry-s2,xc,yc) >= 1 && (rx+s1)*(ry-s2) > rx*ry) {
	rx = rx+s1;
	ry = ry-s2;
	}
	if (check_inner(px,py,rx-s1,ry+s2,xc,yc) >= 1 && (rx-s1)*(ry+s2) > rx*ry) {
	rx = rx-s1;
	ry = ry+s2;
	}

	if (check_inner(px,py,rx+s1,ry,xc+s2,yc) >= 1) { 
	rx = rx+s1;
	xc = xc+s2;
	}
	if (check_inner(px,py,rx+s1,ry,xc-s2,yc) >= 1) {
	rx = rx+s1;
	xc = xc-s2;      
	}
	if (check_inner(px,py,rx,ry+s1,xc,yc+s2) >= 1) {
	ry = ry+s1;
	yc = yc+s2;
	}
    
	if (check_inner(px,py,rx,ry+s1,xc,yc-s2) >= 1) { 
	ry = ry+s1;
	yc = yc-s2;      
	}
	if (check_inner(px,py,rx,ry,xc,yc) < 1) {
	rx = rx-2*step-0.001;
	ry = ry-2*step-0.001;
	if (rx < 0.001)
	rx = 0.001;
	if (ry < 0.001)
	ry = 0.001;

	}
   
	if (rx == rxprev && ry == ryprev && xcprev == xc && ycprev== yc  && check_inner(px,py,rx,ry,xc,yc) >= 1) { 
	vec jxy = get_jumperpoint(px,py,step*5,rx,ry,xc,yc);
	double jx = jxy[0];
	double jy = jxy[1];
      
	xc = xc+jx;
	yc = yc+jy;
	}
	if (rx*ry > maxarea) {
	bestiter = iterat;
	bestxy[0] = xc;
	bestxy[1] = yc;
	besta = rx;
	bestb = ry;
	maxarea = rx*ry;
	}

    }
    return Rcpp::List::create(Named("center") = bestxy,
			      Named("radius.x") = besta,
			      Named("radius.y") = bestb,
			      Named("maxarea") = maxarea);
  
  } catch (std::exception& e) {
    forward_exception_to_r( e );
  } catch (...) {
    ::Rf_error("unknown exception");
  } return R_NilValue; 
} 
