#ifndef POINTS_2_MESH_H_
#define POINTS_2_MESH_H_

#include <RcppArmadillo.h>
#include <cmath>
#include "doozers.h"


using namespace Rcpp;
using namespace arma;

mat updateSearchStruct(mat vb, umat it, uvec clostInd);

double pt_triangle(vec point, vec vbtmp, vec& clost, int& region);

double pt_triplane(vec point, vec vbtmp, vec& clost);

vec pt2mesh(vec point, mat DAT, double& dist, int& faceptr, int& region, int method);

vec getBaryCent(vec point, int fptr, mat vb, umat it);

RcppExport SEXP points2mesh(SEXP ref_,SEXP vb_, SEXP it_, SEXP normals_, SEXP clostInd_, SEXP sign_, SEXP bary_, SEXP method_);

#endif /*POINTS_2_MESH_H_*/
