fx<-function(refmat,M,coefs)
{ q<-dim(M)[1]
  p<-dim(refmat)[1]
  m<-dim(refmat)[2]
  splM<-numeric(q)
    if (m==3)
    {
    for (i in 1:q)
      { Z<-apply((refmat-matrix(M[i,],p,m,byrow=T))^2,1,sum)
        splM[i]<-coefs[p+1]+coefs[p+2]*M[i,1]+coefs[p+3]*M[i,2]+coefs[p+4]*M[i,3]+sum(coefs[1:p]*sqrt(Z))
      }
    }
    else if (m==2)
    {
     for (i in 1:q)
      { Z<-apply((refmat-matrix(M[i,],p,2,byrow=T))^2,1,sum)
	
		Z1<-Z*log(Z)
		Z1[which(is.na(Z1))]<-0
		splM[i]<-coefs[p+1]+coefs[p+2]*M[i,1]+coefs[p+3]*M[i,2]+sum(coefs[1:p]*Z1)
	
      }
    }
    return(splM)
}
