mc.fx<-fx<-function(refmat,M,coefs,time=TRUE)
{ 	
	storage.mode(M)<-"double"
	storage.mode(refmat)<-"double"
	storage.mode(coefs)<-"double"
	q<-dim(M)[1]
  	p<-dim(refmat)[1]
  	m<-dim(refmat)[2]
	
	
	splM<-.Fortran("tpsfx",refmat,p,M,q,refmat[,1],coefs,M[,1])[[7]]
	
	
	
    
    return(splM)
}
