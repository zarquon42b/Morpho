showPC<-function(scores,PC,mshape)
  {
    dims<-dim(mshape)
    if (length(scores) > 1)                                    #pred1<-t(coeff)%*%mod
      predPC<-PC%*%scores
    else
      predPC<-PC*scores
    modell<-mshape+matrix(predPC,dims[1],dims[2])
    return(modell)
}
