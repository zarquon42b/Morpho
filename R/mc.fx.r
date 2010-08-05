mc.fx<-function(refmat,M,coefs,time=TRUE)
{ q<-dim(M)[1]
  p<-dim(refmat)[1]
  m<-dim(refmat)[2]
	aa<-c(1:q)
	aa<-as.list(aa)
	
	
	
	splineFX<-function(x)
	
		{ 
		if (m==3)
			{
			 Z<-apply((refmat-matrix(M[x,],p,m,byrow=T))^2,1,sum)
        		splMi<-coefs[p+1]+coefs[p+2]*M[x,1]+coefs[p+3]*M[x,2]+coefs[p+4]*M[x,3]+sum(coefs[1:p]*sqrt(Z))
      
				
			#Z<-apply((refmat-matrix(x,p,m,byrow=T))^2,1,sum)
        		#splMi<-coefs[p+1]+coefs[p+2]*x[1]+coefs[p+3]*x[2]+coefs[p+4]*x[3]+sum(coefs[1:p]*sqrt(Z))
      			}
		else
			{Z<-apply((refmat-matrix(x,p,2,byrow=T))^2,1,sum)
			Z1<-Z*log(Z)
			Z1[which(is.na(Z1))]<-0
			splMi<-coefs[p+1]+coefs[p+2]*x[1]+coefs[p+3]*x[2]+sum(coefs[1:p]*Z1)
			}
		return(splMi)
		} 
	#if (time)
	#	{t0<-Sys.time()	;test1<-splineFX(M[1,]);t1<-Sys.time();cat(paste("warping....\nestimated time: ",difftime(t1,t0,unit="mins")*q*m)," minutes\n" ,sep="")
		
	splM<-mclapply(aa,splineFX)
	splM<-unlist(splM)
#	splM<-apply(M,1,splineFX)
	
	
	#splM<-numeric(q)
    
    
    return(splM)
}
