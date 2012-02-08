histGroup <- function(data,groups,xlab=dataname,ylab,alpha=0.5, main=paste("Histogram of" , dataname))
  {
   
    dataname <- paste(deparse(substitute(data), 500), collapse="\n")
    
    histo <- hist(data,plot=FALSE)
    

    if(!is.factor(groups))
      {
        groups <- as.factor(groups)
      }
    lev <- levels(groups)
    nlev <- length(lev)
    colo <- rainbow(nlev,alpha=alpha)
    hist(data[groups==lev[1]],breaks=histo$breaks,col=colo[1],main=main,xlab=xlab,ylab=ylab)
    for (i in 2:nlev)
      {
        hist(data[groups==lev[i]],breaks=histo$breaks,col=colo[i],add=T)
      }
  }
