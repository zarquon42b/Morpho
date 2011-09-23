get.obj.size <- function()
  {
    allob <- 0
    wspace <- ls(.GlobalEnv)
    print(wspace)
    for (i in 1:length(wspace))
      {
        allob[i] <- (object.size(get(wspace[i])))
      }
    allob <- data.frame(allob,wspace)
    allob <- allob[order(allob[,1],decreasing=TRUE),]
    return(allob)
  }
