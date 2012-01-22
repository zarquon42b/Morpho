showPC<-function(scores,PC,mshape)
	{dims<-dim(mshape)
	#pred1<-t(coeff)%*%mod
	predPC<-PC%*%scores
	modell<-mshape+matrix(predPC,dims[1],dims[2])
	return(modell)
}
