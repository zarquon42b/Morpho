#include "covPCA.h"
#include "doozers.h"


double covDist(mat &s1, mat &s2) {
  double cdist;
  mat X, EigVec;
  cx_vec eigval, tmp(1);
  bool check = solve(X, s1, s2);
  if (!check)
    return 2.0e30;
  eig_gen(eigval,X);
  eigval = log(eigval);
  tmp = (dot(eigval,eigval));
  cdist = real(tmp(0));
  return(cdist);//squared distance
}

// get pairwise squared distances for multiple groups
mat covDistMulti(mat &data, ivec groups, bool scramble) {
  typedef unsigned int uint;
  uint maxlev = groups.max();
  mat dists(maxlev,maxlev);
  double check;
  List covaList(maxlev);
  // compute covariance matrix for each group
  for (uint i = 0; i < maxlev; ++i) {
    if (!scramble) {
      covaList[i] = cov(data.rows(arma::find(groups == (i+1))));
    } else {//bootstrapping within groups with replacement
      mat tmpdat = data.rows(arma::find(groups == (i+1)));
      uint nrow = tmpdat.n_rows;
      uvec shaker = unique(randi<uvec>(nrow, distr_param(0,nrow-1)));
      //uvec shaker = unique(randi(nrow, 0, nrow-1));
      covaList[i] = cov(tmpdat.rows(shaker));
    }
  }
  // compute pairwise distances between covariance matrices
  dists.zeros();
  for (uint i = 0; i < (maxlev-1); ++i) {
    for (uint j = i+1; j < (maxlev); ++j){
      mat tmp0 = covaList[i];
      mat tmp1 = covaList[j];
      check = covDist(tmp0,tmp1);
      dists(j,i) = check;
    }
  }
  mat checkerr = dists;
  checkerr.reshape(maxlev*maxlev,1);
  colvec checkvec = checkerr.col(0);  
  if (any(checkvec == 2.0e30)) {//check if covDist failed solving
    mat errout(0,0);
    return errout;
  } else {
    dists = dists+dists.t();    
    return dists;
  }
  
}
cube covPCAboot(mat &data, ivec groups, int rounds) {
  typedef unsigned int uint;
  uint maxlev = groups.max();
  mat reference = covDistMulti(data, groups, false);
  List refList = covMDS(reference);
  mat refscores = refList["PCscores"];
  cube alldist(maxlev, maxlev-1, rounds);
  for (int i = 0; i < rounds;){
    mat result = covDistMulti(data, groups, true);
    if (result.n_cols > 0) {
      List tmplist = covMDS(result);
      bool check = tmplist["check"];
      if (check) {
	mat tmpscores = tmplist["PCscores"];
	for (uint j =0; j < tmpscores.n_cols; j++) {
	  if (angcalcArma(tmpscores.col(j),refscores.col(j)) > 1.570796f )
	    tmpscores.col(j) *= -1;
	}
	alldist.slice(i) = tmpscores;
	i++;//only increment if covPCA did not fail
      }
    }
  }
  return alldist;
}

// permutation tests by shuffling group affinities
cube covPCApermute(mat &data, ivec groups, int rounds) {
  typedef unsigned int uint;
  uint maxlev = groups.max();  
  cube alldist(maxlev, maxlev, rounds);
  for (int i = 0; i < rounds;){
    groups = shuffle(groups);
    mat result = covDistMulti(data, groups, false);
    result = sqrt(result);
    if (result.n_cols > 0) {
      alldist.slice(i) = result;
      i++;//only increment if covPCA did not fail
    }
  }
  return alldist;

}
// compute PCscores using MDS approach
List covMDS(mat &dists) {
  unsigned int nlev = dists.n_cols;
  double hf = nlev;
  hf = -1/hf;
  mat H(nlev,nlev);
  H.fill(hf);
  H.diag() += 1;
  mat D = -0.5*(H*dists*H);
  mat eigvec;
  vec eigval;
  mat PCscores;
  bool check = eig_sym(eigval, eigvec, D);
  if (check) {
    uvec useandsort(nlev-1);//sort eigenvectors and values by increasing value
    for (unsigned int i = 0; i < (nlev-1); i++)
      useandsort(i) = nlev-1-i;
    eigval = eigval.elem(useandsort);
    eigvec = eigvec.cols(useandsort);
    
    PCscores = eigvec;
    for (unsigned int i = 0; i < (nlev-1); i++) {
      if (eigval(i) > 0)
	PCscores.col(i) *= sqrt(eigval(i));
      else 
	check = false;
    }
    }
  List out = List::create(Named("eigenvec")=eigvec,
			  Named("eigenval")=eigval,
			  Named("PCscores")=PCscores,
			  Named("check")=check
			  );
  return out;
}

//this is the function exposed to R in covPCA
SEXP covPCAwrap(SEXP data_, SEXP groups_, SEXP scramble_, SEXP rounds_) {
  try {
  //process input
  int scramble = Rcpp::as<int>(scramble_);
  int rounds = as<int>(rounds_);
  mat armaData = as<mat>(data_);
  arma::ivec armaGroups = Rcpp::as<arma::ivec>(groups_);
  
  // get distance matrix
  mat dist = covDistMulti(armaData,armaGroups,false);
  cube bootstrap;
  // run bootstrapping of PCscores.
  if (scramble > 0)
    bootstrap = covPCAboot(armaData,armaGroups,scramble);
  
  List out = covMDS(dist);

  // run permutations and store resulting distances in a 3D-cube
  cube permutest;
  if (rounds > 0)
    permutest = covPCApermute(armaData,armaGroups,rounds);
  //setup List to return to R
  return List::create(Named("dist")=sqrt(dist),
		      Named("Scores")=out,
		      Named("permute")=permutest,
		       Named("bootstrap")=bootstrap
		      );
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
//this is a function exposed to R calculating distance only 
//currently not used  
SEXP covWrap(SEXP s1_, SEXP s2_) {
  mat S1 = as<mat>(s1_);
  mat S2 = as<mat>(s2_);
  double cdist = covDist(S1, S2);
  return wrap(cdist);

}  
