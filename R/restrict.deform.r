restrict <- function(x,model,sd=3,maxVar=95,scale=FALSE)
  {
    dims <- dim(x)
    mshape <- model$mshape
    PCs <- model$PCs
    pc.used <-which(model$Variance[,3] < maxVar)
    sds <-model$Variance[pc.used,1]
    prob <- TRUE
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
                prob=FALSE
                xscore[i] <- sd*sq.sds[i]*signum
              }
          }
      }
    
    restr.x <- matrix(PCs[,pc.used]%*%xscore,dims[1],dims[2])+mshape
    return(list(restr.x=restr.x,prob=prob))
  }

warp.restrict <- function(x,which,tar.lm,model,tol=1e-3,sd=3,maxVar=95,scale=F,recurse=T,uniform=TRUE,iterations=NULL)
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
        tmp <- tmp/c.size(tmp)
        tmp <- rotonto(model$mshape,tmp,scaling=F)$yrot
        tmp.res <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale)
        tmp <-tmp.res$restr.x
        if (tmp.res$prob)
          {
            #if (!sd.i >= sd)
            p <- 0
          }
        tmp <-rotonmat(tmp,tmp[which,],tar.lm,scale=T)
        
        if (recurse)
          {
            tmp.lm <- tmp[which,]
          }
        p <- angle.calc(tmp,tmp.old)$rho
                                        # print(p)
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
    
    clean.out <- tmp
    clean.out[which,] <-tar.lm 
    
    return(list(raw=tmp,clean=clean.out))
  }
