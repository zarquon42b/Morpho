mc.fx<-fx<-function(refmat,M,coefs,time=TRUE)
{ 	
	print(coefs)
	q<-dim(M)[1]
  	p<-dim(refmat)[1]
  	m<-dim(refmat)[2]
	if(m==3)
		{M1<-cbind(1,M)
		coefs<-t(coefs)
		storage.mode(M)<-"double"
		storage.mode(refmat)<-"double"
		storage.mode(coefs)<-"double"
		splM<-.Fortran("tpsfxo",refmat,p,M,q,M1,refmat[,1],coefs,M)[[8]]
		}
	
	else if (m==2)
    {splM<-M

     for (i in 1:q)
      { ### calculate nonaffine vector
		
		Z<-apply((refmat-matrix(M[i,],p,2,byrow=T))^2,1,sum)
		Z1<-Z*log(Z)
		Z1[which(is.na(Z1))]<-0
		M1<-cbind(1,M)
		print(M1)
		splM[i,]<-t(coefs[(p+1):(p+3),])%*%M1[i,]+t(coefs[1:p,])%*%Z1
	
      }
	}


	
    
    return(splM)
}
