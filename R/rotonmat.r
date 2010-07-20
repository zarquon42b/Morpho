rotonmat<-function(X,refmat,tarmat,scale=TRUE)
{	
	ro<-rotonto(tarmat,refmat,scaling=scale,signref=F)
	Xrot<-t(apply(X,1,function(x){x-ro$transy}))%*%ro$gamm
	return(Xrot)
}
