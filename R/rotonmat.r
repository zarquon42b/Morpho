rotonmat<-function(X,refmat,tarmat,scale=TRUE)
{	
	ro<-rotonto(tarmat,refmat,scale=scale,signref=F)
	
	Xrot<-t(apply(X,1,function(x){x-ro$transy}))%*%ro$gamm
	if (scale)
		{sf<-cSize(refmat)/cSize(ro$yrot)
		Xrot<-Xrot/sf
		}
	Xrot<-t(apply(Xrot,1,function(x){x+ro$trans}))
	return(Xrot)
}
