procAOV <- function(symproc,indnames=NULL)
  {
    
    if (class(symproc) != "symproc")
      {stop("input is not of class symproc")
     }
    m <- dim(symproc$rotated)[2]
    k <- dim(symproc$rotated)[1]
    n <- dim(symproc$rotated)[3]

    if (is.null(indnames))
      indnames <- as.factor(rownames(symproc$PCscore_sym))
    else
      indnames <- as.factor(indnames)
    indlev <- levels(indnames)
    nlev <- length(indlev)
    alist <- list()
    
    for (i in 1:length(indlev))
      {
       alist[[i]] <- grep(indlev[i],indnames)
      }

    checkinds <- lapply(alist,function(x){length(x) == length(alist[[1]])})
    if (prod(as.integer(unlist(checkinds))) == 0)
      stop("same number of digitizations is needed for all specimen")
    r <- length(alist[[1]])
    pl <- dim(symproc$pairedLM)[1]
    sl <- k-2*pl

    df.ind <- (nlev-1)*(3*pl+2*sl-4)
    df.side <- 3*pl+sl-3
    df.indxside <- (nlev-1)*(3*pl+sl-3)
    df.res <- (r-1)*nlev*(3*(2*pl+sl)-7)
    side <- sum(symproc$asymmean^2)*n
     

    asymfun <- function(x)
      {
        x <- sum(apply(symproc$Asymtan[x,],2,mean)^2)*length(x)
        return(x)
      }
    symfun <- function(x)
      {
        x <- sum(apply(symproc$Symtan[x,],2,mean)^2)*length(x)
        return(x)
      }
    if (length(alist[[1]]) > 1)
        {
          
          indxside <- sum(unlist(lapply(alist,asymfun)))
          ind <- sum(unlist(lapply(alist,symfun)))
        }
        else
        {indxside <- sum(symproc$Asymtan^2)
         ind <- sum(symproc$Symtan^2)
        }
    allsq <- sum(symproc$Asymtan[,]^2+symproc$Symtan[,]^2)+side
    res <- allsq-(ind+indxside+side)
    
        
    outss <- c(ind,side,indxside,res)
    outdf <- as.integer(c(df.ind,df.side,df.indxside,df.res))
    exVar <- outss/allsq
    outms <- outss/outdf
    out <- data.frame(outss,outms,exVar,outdf)
    F.values <- c(outms[1]/outms[3],outms[2]/outms[3],outms[3]/outms[4],NA)
    sig <- c(pf(F.values[1],outdf[1],outdf[3],lower.tail=F),pf(F.values[2],outdf[2],outdf[3],lower.tail=F),pf(F.values[3],outdf[3],outdf[4],lower.tail=F),NA);

   
    out <- data.frame(out,F.values,sig)
    rownames(out) <- c("ind","side","ind.x.side","error")
    
    colnames(out) <- c("SS","MS","exVar","df","F","p-value")
    return(out)
  }
