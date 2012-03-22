name2factor <- function(x,sep="_",which,collapse=sep,dif=NULL)
  {
    if (length(dim(x))==3)
      {
        names <- dimnames(x)[[3]]
      }
    else if (length(dim(x))==2)
      {
        names <- dimnames(x)[[1]]
      }
    else
      {
        names <- names(x)
      }
    fac <- strsplit(names,split=sep)
    
    for (i in 1:length(fac))
      {
                 
        fac[[i]] <- paste(fac[[i]][which],collapse=collapse)
        
      }
    fac <- as.factor(unlist(fac))
    return(fac)
  }
name2num <- function(x,sep="_",which,collapse=sep,dif=TRUE)
  {
    if (length(dim(x))==3)
      {
        names <- dimnames(x)[[3]]
      }
    else if (length(dim(x))==2)
      {
        names <- dimnames(x)[[1]]
      }
    else
      {
        names <- names(x)
      }
    fac <- strsplit(names,split=sep)
    
    for (i in 1:length(fac))
      {
        if (dif)
          {
            fac[[i]] <- diff(as.numeric(fac[[i]][which]))
          }
        else
          {
            fac[[i]] <- as.numeric(fac[[i]][which])
          }
      }
    fac <-(unlist(fac))
    return(fac)
  }
