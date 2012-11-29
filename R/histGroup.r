histGroup <- function(data,groups, main=paste("Histogram of" , dataname),xlab=dataname,ylab,alpha=0.5,breaks="Sturges")
  {
   
    dataname <- paste(deparse(substitute(data), 500), collapse="\n")
    
    histo <- hist(data,plot=FALSE,breaks=breaks)
    

    if(!is.factor(groups))
      {
        groups <- as.factor(groups)
      }
    lev <- levels(groups)
    nlev <- length(lev)
    colo <- rainbow(nlev,alpha=alpha)
    testrun <- 0
    for( i in 1:nlev)
     { testrun[i] <-  max(hist(data[groups==lev[i]],breaks=histo$breaks,plot=F)$counts)
     }
   ylim <- max(testrun)
    ylim <- ylim+0.15*ylim
    hist(data[groups==lev[1]],breaks=histo$breaks,col=colo[1],main=main,xlab=xlab,ylab=ylab,ylim=c(0,ylim))
    for (i in 2:nlev)
      {
        hist(data[groups==lev[i]],breaks=histo$breaks,col=colo[i],add=T)
      }
  }
