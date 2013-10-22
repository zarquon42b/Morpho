write.pts <- function(x, filename=dataname)
{
    dataname <- deparse(substitute(x))
    filename <- paste(filename,".pts",sep="")
    k <- dim(x)[1]
    m <- dim(x)[2]
    x <- as.matrix(x)
    a0 <- paste("S",000,c(0:(k-1)),sep="")
    all.frame <- cbind(a0,x)
    cat("Version 1.0\n",file=filename)
    cat(paste(k,"\n",sep=""),file=filename,append=T)
    write(t(all.frame),file=filename,sep=" ",ncolumns =  m+1,append=TRUE)
}
