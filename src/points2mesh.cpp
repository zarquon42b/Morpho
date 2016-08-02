#include "points2mesh.h"

// update search structures only taking into acount probable face candidates
mat updateSearchStruct(mat vb, umat it, uvec clostInd) {
  int nit = it.n_cols;
  int used = clostInd.size();
  mat DAT(13,nit);
  DAT.zeros();
  for (int i=0; i < used; ++i) {
    int pos = clostInd(i);
    uvec ptr = linspace<uvec>(0,2,3);
    uvec ui(1); ui.fill(pos);
    uvec ittmp = it.col(pos);
    mat vbtmp = vb.cols(ittmp);
    DAT(ptr,ui) = vbtmp.col(0);//B1
    DAT(ptr+3,ui) = vbtmp.col(1)-vbtmp.col(0);//e0
    DAT(ptr+6,ui) =  vbtmp.col(2)-vbtmp.col(0);//e1
    DAT(9,pos) = dot(DAT(ptr+3,ui),DAT(ptr+3,ui));//a
    DAT(10,pos) = dot(DAT(ptr+3,ui),DAT(ptr+6,ui));//b
    DAT(11,pos) = dot(DAT(ptr+6,ui),DAT(ptr+6,ui));//c
    DAT(12,pos) = std::abs(DAT(9,pos)*DAT(11,pos) - pow(DAT(10,pos),2));
  }
  return(DAT);
}
// search clostest point on triangle
double pt_triangle(vec point, vec vbtmp, vec& clost, int& region) {
  double sqdist;
  uvec ptr = linspace<uvec>(0,2,3);
  vec B = vbtmp(ptr);
  vec dv = B - point;
  vec e0 = vbtmp(ptr+3);
  vec e1 = vbtmp(ptr+6);
  double a00 = vbtmp(9);
  double a01 = vbtmp(10);
  double a11 = vbtmp(11);
  double b0 = dot(e0,dv);
  double b1 = dot(e1,dv);
  double c = dot(dv,dv);
  double det =  vbtmp(12);
  double s = a01*b1- a11*b0;
  double t = a01*b0 - a00*b1;
  double invDet, tmp0, tmp1, numer ,denom;
  if (s+t <= det) {
    if (s < 0) {
      if (t < 0) {
	//region 4 begin
	region = 4;
	if (b0 < 0) {
	  t = 0;
	  if ( -b0 >= a00) {
	    s = 1;
	    sqdist = a00 + 2*b0 + c;
	  } else {
	    s = -b0/a00;
	    sqdist = b0*s + c;
	  }
	} else {
	  s = 0;
	  if ( b1 >= 0 ) {
	    t=0;
	    sqdist = c;
	  } else if (-b1 >= a11) {
	    t = 1;
	    sqdist = a11 + 2*b1 + c;
	  } else {
	    t = -b1/a11;
	    sqdist = b1*t +c;
	  }
	}
	//region 4 end
      } else {
	//region 3 begin
	region = 3;
	s = 0;
	if (b1 >=0) {
	  t=0;
	  sqdist = c;
	} else if (-b1 >= a11) {
	  t = 1;
	  sqdist = a11+2*b1+c;
	} else {
	  t = -b1/a11;
	  sqdist = b1*t+c;
	}
	//region3 end
      }
    } else if (t < 0) {
      //region 5 begin
      region = 5;
      t = 0;
      if (b0 >= 0) {
	s = 0;
	sqdist = c;
      } else if (-b0 >=a00) {
	s = 1;
	sqdist = a00 + 2*b0 + c;
      } else {
	s = -b0/a00;
	sqdist = b0*s + c;
      }
      //region 5 end
    } else {      
      //region 0 begin
      region = 0;
      invDet = 1/det;
      s = s*invDet;
      t=  t*invDet;
      sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
      // region 0 end
    }
  } else {
    if (s < 0) {
      //region 2 begin
      region = 2;
      tmp0 = a01 + b0;
      tmp1 = a11 + b1;
      if (tmp1 > tmp0) {
	numer = tmp1 - tmp0;
	denom = a00 - 2*a01 + a11;
	if (numer >= denom) {
	  s = 1;
	  t = 0;
	  sqdist = a00 + 2*b0 + c;
	} else {               
	  s = numer/denom;
	  t = 1 - s;
	  sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t +2*b1) + c;
	}
      } else { 
	s = 0;
	if (tmp1 <= 0) {
	  t = 1;
	  sqdist = a11+2*b1+c;
	} else if (b1 >= 0) {
	  t = 0;
	  sqdist = c;
	} else {
	  t = -b1/a11;
	  sqdist = b1*t + c;
	}
	//region 2 end
      }
    } else if (t < 0) {
      //bin/region 6 begin
      region = 6;
      tmp0 = a01 + b1;
      tmp1 = a00 + b0;
      if (tmp1 > tmp0) {
	numer = tmp1 - tmp0;
	denom = a00 - 2*a01 + a11;
	if (numer >= denom) {                             
	  t = 1;
	  s = 0;
	  sqdist = a11 + 2*b1 + c;
	} else {
	  t = numer/denom;
	  s = 1 - t;
	  sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
	}
      } else {
	t = 0;
	if (tmp1 <= 0) {
	  s = 1;
	  sqdist = a00 + 2*b0 + c;
	} else if ( b0 >=0) {
	  s = 0;
	  sqdist = c;
	} else {
	  s = -b0/a00;
	  sqdist = b0*s +c;
	}
      }
    } else {
      numer = a11 + b1 - a01 - b0;
      if (numer <= 0) {
	s = 0;
	t = 1; 
	sqdist = a11+ 2*b1 + c;
      } else {
	denom = a00 - 2*a01 + a11;;
	if (numer >= denom) {
	  s =1;
	  t = 0;
	  sqdist = a00 + 2*b0 + c;
	} else {
	  s = numer/denom;
	  t = 1 - s;
	  sqdist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
	}
      }
    }
  }
  clost = B + s*e0 + t*e1;
  return sqdist;
}

// alternate method including face orientation as by Moshfeghi
double pt_triplane(vec point, vec vbtmp, vec& clost) {
  double sqdist;
  uvec ptr = linspace<uvec>(0,2,3);
  vec B = vbtmp(ptr);
  vec dv = B - point;
  vec e0 = vbtmp(ptr+3);
  vec e1 = vbtmp(ptr+6);
  vec normal(3); normal.zeros();
  crosspArma(e0,e1,normal);
  double normlen = norm(normal,2);
  if (normlen > 0)
    normal /=normlen;
  vec diffvec = point-B;
  double difflen = norm(diffvec,2);
  if (difflen > 0)
    diffvec /= difflen;
  double alpha = dot(diffvec,normal);
  double p0p0dist = sqrt(dot(B - point,B-point))*alpha;
  vec p0p0 = -p0p0dist*normal;
  clost = point + p0p0;
  sqdist = p0p0dist*p0p0dist;
  return sqdist;

}


vec pt2mesh(vec point, mat DAT, double& dist, int& faceptr, int& region, int method) {
  int ndat = DAT.n_cols;
  vec closttmp(3); closttmp.fill(9999);
  vec clost(3);
  vec checkclost(3);
  vec vbtmp(13);
  double dist_old = 1e10;
  int regiontmp = 0;
  double sqdist = 0;
  bool meth = false;
  for (int i=0; i < ndat; ++i) {
    vbtmp = DAT.col(i);
    sqdist = pt_triangle(point, vbtmp, closttmp, regiontmp);
    if (method == 1 && regiontmp != 0 ) {
      meth = true;
      sqdist = pt_triplane(point,vbtmp,checkclost);
      sqdist = sqrt(sqdist) + sqrt(dot(checkclost-closttmp,checkclost-closttmp));
      sqdist = sqdist*sqdist;
    }
    if (sqdist < dist_old) {
      dist_old = sqdist;
      clost = closttmp;
      region = regiontmp;
      faceptr = i;
      if (meth) // get correct distance
      	dist = norm(point-clost,2);
      else 
	dist = sqrt(std::abs(sqdist));
    }
  }
  return clost;
}
// calculate barycentric coordinates of a point, given a face pointer and a mesh.
vec getBaryCent(vec point, int fptr, mat vb, umat it) {
  vec barycoord(3); barycoord.zeros();
  vec v0 = vb.col(it(1,fptr))-vb.col(it(0,fptr));
  vec v1 = vb.col(it(2,fptr))-vb.col(it(0,fptr));
  vec v2 = point - vb.col(it(0,fptr));
  double d00 = dot(v0,v0);
  double d01 = dot(v0,v1);
  double d11 = dot(v1,v1);
  double d20 = dot(v2,v0);
  double d21 = dot(v2,v1);
  double denom = d00*d11-d01*d01;
  barycoord(1) = (d11 * d20 - d01 * d21) / denom;
  barycoord(2) =  (d00 * d21 - d01 * d20) / denom;
  barycoord(0) = 1 - barycoord(1) - barycoord(2);
  return barycoord;
}
// main function to handle in and output
SEXP points2mesh(SEXP ref_,SEXP vb_, SEXP it_, SEXP normals_, SEXP clostInd_, SEXP sign_, SEXP bary_, SEXP method_) {
  try {
    mat ref = as<mat>(ref_);//reference 
    mat vb = as<mat>(vb_);//target vertices
    mat normals = as<mat>(normals_);//target normals
    umat clostIndU = as<arma::umat>(clostInd_);//face indices to search on
    int nref = ref.n_cols;
    bool sign = as<bool>(sign_);
    bool bary = as<bool>(bary_);
    int method = as<int>(method_);
    umat itU = as<arma::umat>(it_);// convert to unsigned
    // check which faces are acutally searched
    uvec uniclost = unique(clostIndU);
    // calculate edges and stuff needed for point search
    mat DAT = updateSearchStruct(vb,itU,uniclost);
    mat closeMat = ref;
    mat outnormals = ref;
    mat barycoords = ref;
  
    ivec region(nref), faceptr(nref); region.zeros();faceptr.zeros();
    vec dists(nref); dists.zeros();
    for (int i=0; i < nref; ++i) {
      mat tmpdat = DAT.cols(clostIndU.col(i));//select appropriate subset of DAT
      vec tmpvec(3), weight(3); 
      tmpvec.zeros();weight.zeros();
      int faceptrtmp=0;
      int regioni = region(i);
      closeMat.col(i) = pt2mesh(ref.col(i), tmpdat, dists(i), faceptrtmp, regioni,method);
      faceptr(i)=clostIndU(faceptrtmp,i);
      //get normal weights
      for (int j =0; j < 2; ++j) {
	vec tmpdiff = closeMat.col(i)-vb.col(itU(j,faceptr(i)));
	weight(j) = sqrt(dot(tmpdiff,tmpdiff));
	if (weight(j) == 0.0) {
	  weight(j) = 1e12;
	} else {
	  weight(j) = 1/weight(j);
	}
      }
      //get weighted normals
      vec tmpnormal(3); tmpnormal.zeros();

      for (int j = 0; j < 3; j++) {
	tmpnormal += weight(j)*normals.col(itU(j,faceptr(i)));
      }
      double normlen = norm(tmpnormal,2);
      if (normlen > 0) {
	tmpnormal /= normlen;
      }
      outnormals.col(i) = tmpnormal;
      // sign distances
      if (sign) {
	vec tmpdiff = closeMat.col(i) - ref.col(i);
	double signo = dot(tmpdiff,tmpnormal);
	if (signo < 0) 
	  dists(i) *= -1;
      }
      // get barycentric coords
      if (bary)
	barycoords.col(i) = getBaryCent(closeMat.col(i), faceptr(i), vb, itU);
    }  
    return Rcpp::List::create(Named("clost") = closeMat,
			      Named("dists") = dists,
			      Named("faceptr") = faceptr,
			      Named("normals") = outnormals,
			      Named("barycoords") = barycoords
			      );
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
  
  
}



  
