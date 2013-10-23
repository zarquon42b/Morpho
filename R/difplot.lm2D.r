.difplotLM2D <- function(refshape,targetshape,color=4,lwd=1,lcol=2,main=main, text=TRUE)
{   
    A <- refshape
    k <- dim(A)[1]
    m <- dim(A)[2]
    sds <- 0
    plot(refshape,  col = color,main=main,asp=1,axes=FALSE,xlab="",ylab="")
    for (j in 1:k)
        lines(rbind(refshape[j,],targetshape[j,]),col=lcol,lwd=lwd)
    if (text)
        text(refshape,labels=paste("",c(1:k),sep=""),col=lcol,pos=2) 
}
