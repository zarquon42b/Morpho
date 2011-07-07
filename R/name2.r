name2factor <- function(x,sep="_",which,collapse=sep,dif=NULL)
  {
    names <- dimnames(x)[[3]]
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
    names <- dimnames(x)[[3]]
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
