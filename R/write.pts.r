#' exports a matrix containing landmarks into .pts format
#' 
#' exports a matrix containing landmarks into .pts format that can be read by
#' IDAV Landmark.
#' 
#' you can import the information into the program landmarks available at
#' http://graphics.idav.ucdavis.edu/research/EvoMorph
#' 
#' @param x k x m matrix containing landmark configuration
#' @param filename character: Path/name of the requested output - extension
#' will be added atuomatically. If not specified, the file will be named as the
#' exported object.
#' @author Stefan Schlager
#' @seealso \code{\link{read.pts}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data(nose)
#' write.pts(shortnose.lm, filename="shortnose")
#' 
#' @export
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
