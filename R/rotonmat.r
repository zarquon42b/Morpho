rotonmat<-function(X,refmat,tarmat,scale=TRUE,scaling=FALSE,reflection=FALSE)
{	
	ro<-rotonto(tarmat,refmat,scale=scale,signref=F,reflection=reflection)
	
	Xrot<-t(apply(X,1,function(x){x-ro$transy}))%*%ro$gamm
	if (scale)
		{sf<-cSize(refmat)/cSize(ro$yrot)
		Xrot<-Xrot/sf
		}
	Xrot<-t(apply(Xrot,1,function(x){x+ro$trans}))
	return(Xrot)
}
