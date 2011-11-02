restrict <- function(x,model,sd=3,maxVar=95,scale=FALSE)
  {
    dims <- dim(x)
    mshape <- model$mshape
    PCs <- model$PCs
    pc.used <-which(model$Variance[,3] < maxVar)
    sds <-model$Variance[pc.used,1]
    
    sdl <- length(pc.used)
    
    xscore <- t(PCs[,pc.used])%*%as.vector(x-mshape)
    if (scale)
      {
        Mt <- qchisq((pnorm(sd)),df=sdl)
        prob <- sum(xscore^2/sds)
        if (prob > Mt )
        {
          sca <- Mt/prob
          xscore <- xscore*sqrt(sca)
        }
      }
    else
      {
        sq.sds <- sqrt(sds)
        for (i in 1:length(xscore))
          {
            signum <- sign(xscore[i])
            if (abs(xscore[i]) > (sd*sq.sds[i]))
              {
                xscore[i] <- sd*sq.sds[i]*signum
              }
          }
      }
    
    restr.x <- matrix(PCs[,pc.used]%*%xscore,dims[1],dims[2])+mshape
    return(restr.x)
  }

warp.restrict <- function(x,which,tar.lm,model,tol=1e-5,sd=3,maxVar=95,scale=F,recurse=T,uniform=TRUE,iterations=NULL)
  {
    x.lm <- x[which,]
    sd.i <- sd
    tmp <- x
    tmp.lm <- x.lm
    p <- 1e10
    count <- 0
    while(p > tol)
      { count <- count+1
        p.old <- p
        tmp.old <- tmp     
        tmp <- tps3d(tmp,tmp.lm,tar.lm)
        tmp <- rotonto(model$mshape,tmp)$yrot        
        tmp <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale)
        if (recurse)
          {
            tmp.lm <- tmp[which,]
          }
        p <- angle.calc(tmp,tmp.old)$rho
       print(p)
        if (uniform)
          {
            if (p > p.old)
              {
                p <- 0
                tmp <- tmp.old
              }                  
          }
        if (!is.null(iterations))
          if (count == iterations)
            {
              p <- 0
            }
      }
    
    return(tmp)
  }
