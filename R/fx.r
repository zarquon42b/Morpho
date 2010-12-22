mc.fx<-fx<-function(refmat,M,coefs,time=TRUE)
{ 	
	
	q<-dim(M)[1]
  	p<-dim(refmat)[1]
  	m<-dim(refmat)[2]
	if(m==3)
		{storage.mode(M)<-"double"
		storage.mode(refmat)<-"double"
		storage.mode(coefs)<-"double"
		splM<-.Fortran("tpsfx",refmat,p,M,q,refmat[,1],coefs,M[,1])[[7]]
		}
	
	else if (m==2)
    {splM<-numeric(q)

     for (i in 1:q)
      { Z<-apply((refmat-matrix(M[i,],p,2,byrow=T))^2,1,sum)
	
		Z1<-Z*log(Z)
		Z1[which(is.na(Z1))]<-0
		splM[i]<-coefs[p+1]+coefs[p+2]*M[i,1]+coefs[p+3]*M[i,2]+sum(coefs[1:p]*Z1)
	
      }
	}


	
    
    return(splM)
}
