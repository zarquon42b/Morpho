get.obj.size <- function(units=2)
  {
    units <- units[1]
    
    allob <- 0
    wspace <- ls(.GlobalEnv)
    #print(wspace)
    for (i in 1:length(wspace))
      {
        allob[i] <- (object.size(get(wspace[i])))
      }
    allob <- round(allob/1024^units,digits=2)
    allob <- data.frame(allob,wspace)
    allob <- allob[order(allob[,1],decreasing=TRUE),]
    return(allob)    
  }
