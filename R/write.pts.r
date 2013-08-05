write.pts <-  function(X,filename)
{   filename <- paste(filename,".pts",sep="")
    k <- dim(X)[1]
    m <- dim(X)[2]
    X <- as.matrix(X)
    a0 <- paste("S",000,c(0:(k-1)),sep="")
    all.frame <-  cbind(a0,X)
    cat("Version 1.0\n",file=filename)
    cat(paste(k,"\n",sep=""),file=filename,append=T)
    write(t(all.frame),file=filename,sep=" ",ncolumns =  m+1,append=TRUE)

    }
