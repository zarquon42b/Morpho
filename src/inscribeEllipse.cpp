#include <RcppArmadillo.h>
#include "inscribeEllipse.h"

using namespace Rcpp;
using namespace std;
using namespace arma;

double check_inner(std::vector<float> poly_x, std::vector<float> poly_y, float a, float b, float x0, float y0) {
  
  float minrr = 5;
  for (int i=0; i < poly_x.size();i++) {
    float x = poly_x[i];
    float y = poly_y[i];
    float rr = pow((x-x0)/a,2) + pow((y-y0)/b,2);
    if (rr < minrr)
      minrr = rr;
  }
  return minrr;
}

vec get_jumperpoint(std::vector<float> poly_x, std::vector<float> poly_y, float sizeup,float a, float b, float x0,float y0) {
  std::vector<double> badx ,bady;
    
  for (int i=0; i < poly_x.size();i++) {
    float x = poly_x[i];
    float y = poly_y[i];
    float rr = pow((x-x0)/(a+sizeup),2)+pow(((y-y0)/(b+sizeup)),2);
    if (rr < 1) {
      badx.push_back(x);
      bady.push_back(y);
    }
  }
  vec badxA(badx);
  vec badyA(bady);
  float xmean = mean(badxA);
  float ymean = mean(badyA);
  vec out(2);
  out[0] = x0-xmean;
  out[1] = y0-ymean;
 
  
  return out;
};




RcppExport SEXP inscribeEllipseCpp(SEXP polyMat_, SEXP step_, SEXP iters_,SEXP init_point_) {
  
  try{
    mat polyMat = as<mat>(polyMat_);
    std::vector<double> px_old =  conv_to<std::vector<double>>::from(polyMat.col(0));
    std::vector<double> py_old = conv_to<std::vector<double>>::from(polyMat.col(1));
    float step = as<float>(step_);
    int iters = as<int>(iters_);
    std::vector<double> init_point = as<std::vector<double>>(init_point_);
    std::vector<double> bestxy = init_point;
    float init_radius = step;
  
    std::vector<float> px, py;
    for (int i = 0; i < (px_old.size()-1); i++) {
      float dx = px_old[i+1]-px_old[i];
      float dy = py_old[i+1]-py_old[i];
      float len1 = pow((dx*dx+dy*dy),0.5);
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
    float minrr = 6;
    float rx = init_radius;
    float besta = init_radius;
    float bestb = init_radius;
    float ry = init_radius;
    float bestx = init_point[0];
    float xc = init_point[0];
    float besty = init_point[1];
    float yc = init_point[1];
    float maxarea = step*step;
    float bestiter = 0;
   
 
    for (int iterat = 0; iterat < iters; iterat++) {
   
      float s1 = step;
      float s3 = step;
      float s2 = step;
      float rxprev = rx;
      float ryprev = ry;
      float xcprev = xc;
      float ycprev = yc;
     

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

RcppExport SEXP inscribeEllipseRotCpp(SEXP polyList_, SEXP step_, SEXP iters_,SEXP init_point_) {
  List polyList(polyList_);
  int n = polyList.size();
  List out;
  float maxarea = 0;
  for (int i = 0; i < n; i++) {
    List temp = inscribeEllipseCpp(polyList[i], step_, iters_,init_point_);
    float tmparea = as<float>(temp["maxarea"]);
    
    if (tmparea > maxarea) {
      maxarea = tmparea;
      out = temp;
      out["bestiter"] = i+1;
    }
  }
  return out;
  
}
  
