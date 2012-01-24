restrict <- function(x,model,sd=3,maxVar=95,scale=FALSE,nPC=NULL)
  {
    dims <- dim(x)
    mshape <- model$mshape
    PCs <- model$PCs
    if (is.null(nPC))
      {
        pc.used <-which(model$Variance[,3] < maxVar)
      }
    else
      {
        pc.used <- 1:nPC
      }
    sds <-model$Variance[pc.used,1]
    prob <- TRUE
    sdl <- length(pc.used)
    
    xscore <- t(PCs[,pc.used])%*%as.vector(x-mshape)

    if (scale)
      {
        Mt <- qchisq(1-2*pnorm(sd,lower.tail=F),df=sdl)
        probs <- sum(xscore^2/sds)
      #  print(Mt)
      #  print(probs)
        if (probs > Mt )
        {
          prob=FALSE
          sca <- Mt/probs
          
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
                #print(i)
                prob=FALSE
                xscore[i] <- sd*sq.sds[i]*signum
              }
          }
      }
    
    restr.x <- matrix(PCs[,pc.used]%*%xscore,dims[1],dims[2])+mshape
    return(list(restr.x=restr.x,prob=prob))
  }

warp.restrict <- function(x,which,tar.lm,model,tol=1e-5,sd=3,maxVar=95,scale=F,recurse=T,uniform=TRUE,iterations=NULL,nPC=NULL,stop.prob=TRUE)
  {
    if (is.null(nPC))
      {
        pc.used <-which(model$Variance[,3] < maxVar)
        print(paste("First ",max(pc.used)," PCs used"))
      }
      
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
        tmp <- tmp/cSize(tmp)
        tmp <- rotonto(model$mshape,tmp,scale=T)$yrot
        tmp.res <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale,nPC=nPC)
        tmp <-tmp.res$restr.x
        
        tmp <-rotonmat(tmp,tmp[which,],tar.lm,scale=T)
        
        if (recurse)
          {
            tmp.lm <- tmp[which,]
          }
        if (tmp.res$prob && stop.prob)
          {
            cat(paste("probable shape within the boundaries of",sd, "sd reached\n"))
            #if (!sd.i >= sd)
            p <- 0
          }
        else
          {
            p <- angle.calc(tmp,tmp.old)$rho
          }                         # print(p)
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
    if (tmp.res$prob)
      { cat(paste("a probable shape was reached after",count,"iterations\n"))
        
      }
    else
      { cat("progress terminated without reaching probability\n")
      }
    tmp.res <- rotonmat(tmp.res$restr.x,tmp.res$restr.x[which,],tar.lm,scale=T)
    clean.out <- tmp
    clean.out[which,] <-tar.lm 
    
    return(list(raw=tmp,clean=clean.out,tmp=tmp.res))
  }
