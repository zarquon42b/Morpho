restrict <- function(x,model,sd=3,maxVar=95,scale=FALSE,nPC=NULL,probab=FALSE,reference=NULL)
  {
    dims <- dim(x)
    mshape <- model$mshape
    PCs <- model$PCs
    restr.x <- NULL
    if (is.null(nPC)) ### select first # of PCs below threshold of total variance
      {
        pc.used <-which(model$Variance[,3] < maxVar)
      }
    else
      {
        pc.used <- 1:nPC ## use predefined # of PCs
      }
    sds <-model$Variance[pc.used,1]
    prob <- TRUE
    sdl <- length(pc.used)
    xtmp <- x-mshape
    if (!is.null(reference))
      {       
        xtmp[reference,] <- 0 ## set reference=mshape
      }
    
    xscore <- t(PCs[,pc.used])%*%as.vector(xtmp)

    if (scale) ### use chisquare distribution of mahalanobis distance
      {
        Mt <- qchisq(1-2*pnorm(sd,lower.tail=F),df=sdl)
        probs <- sum(xscore^2/sds)
        if (probs > Mt )
        {
          prob=FALSE
          sca <- Mt/probs
          xscore <- xscore*sqrt(sca)
        }
      }
    else ### use probability hypercuboid
      {
        sq.sds <- sqrt(sds)        
        for (i in 1:length(xscore)) ## check if PCscores are below threshold
          {
            signum <- sign(xscore[i])
            if (abs(xscore[i]) > (sd*sq.sds[i]))
              {
                prob=FALSE
                xscore[i] <- sd*sq.sds[i]*signum
              }
          }
      }
    if (!probab)
      {
        restr.x <- matrix(PCs[,pc.used]%*%xscore,dims[1],dims[2])+mshape
        if (!is.null(reference))
          {
            restr.x[reference,] <- x[reference,]
          }
      }
    return(list(restr.x=restr.x,prob=prob))
  }

warpRestrict <- function(x,which,tar.lm,model,tol=1e-5,sd=3,maxVar=95,scale=F,recurse=T,uniform=TRUE,iterations=NULL,nPC=NULL,stop.prob=TRUE,spline=TRUE,useReference=FALSE)
  {

    reference <- NULL
    if (useReference)### "which" will be ignored when calculating probability
      {
        reference <- which
      }
    tmp.res <- list()
    tmp.res$prob <-  FALSE
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
    tmp.orig <- tmp
    tmp.old <- tmp     
    tmp <- tps3d(tmp,tmp.lm,tar.lm)## warp onto target
    tmp <- rotonto(model$mshape,tmp,scale=T)$yrot ### register in database space
    
    while(p > tol)
      {
        
        cat(paste("running iteration",count,"\n"))
        p.old <- p
        tmp.old <- tmp
        prob <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale,nPC=nPC,probab=T)$prob
        if (!prob) ### not yet probable
          {
            ## restrict to boundaries
           
            tmp <- restrict(tmp,model=model,sd=sd.i,maxVar=maxVar,scale=scale,nPC=nPC,probab=F,reference=reference)$restr.x
           # print(dim(tmp))
             tmp.orig <- tmp
            if (spline) ###use restricted data and warp it onto target
              {
                tmp.lm <- tmp[which,] 
                tmp <- tps3d(tmp,tmp.lm,tar.lm)
                tmp <- rotonto(model$mshape,tmp,scale=T)$yrot ### register in database space
              }
            else ### replace restricted reference lm with original ones
              {
                tmp[which,] <- tmp.orig[which,]
                tmp <- rotonto(model$mshape,tmp,scale=T)$yrot ### register in database space
              }
          }
        
        if (prob && stop.prob)
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
        count <- count+1
      }
    if (prob)
      {
        cat(paste("a probable shape was reached after",count-1,"iterations\n"))
      }
    else
      {
        cat("progress terminated without reaching probability\n")
      }
    tmp.orig <- rotonmat(tmp.orig,tmp.orig[which,],tar.lm,scale=T)
    tmp <- rotonmat(tmp,tmp[which,],tar.lm,scale=T)
    clean.out <- tmp
    
    return(list(raw=tmp.orig,clean=clean.out,tmp=tmp.res))
  }
