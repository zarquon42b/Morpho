pred.PC<-function(coeff,mod,PC,mshape)
{	dims<-dim(mshape)
	pred1<-t(coeff)%*%mod
	predPC<-PC%*%pred1
	modell<-mshape+matrix(predPC,dims[1],dims[2])
	return(modell)
}
