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
#' @param rownames provide an optional character vector with rownames
#' @author Stefan Schlager
#' @seealso \code{\link{read.pts}}
#' 
#' @examples
#' 
#' data(nose)
#' write.pts(shortnose.lm, filename="shortnose")
#' 
#' @export
write.pts <- function(x, filename=dataname,rownames=NULL)
{
    dataname <- deparse(substitute(x))
    filename <- paste(filename,".pts",sep="")
    k <- dim(x)[1]
    m <- dim(x)[2]
    x <- as.matrix(x)
    if (is.null(rownames))
        a0 <- paste0("S",sprintf("%04d", 0:(k-1)))
    else
        a0 <- rownames
    all.frame <- cbind(a0,x)
    all.frame <- format(all.frame,trim=TRUE)
    cat("Version 1.0\n",file=filename)
    cat(paste(k,"\n",sep=""),file=filename,append=T)
    write.table(all.frame,append = T,file = filename,col.names = FALSE,quote=FALSE,row.names = FALSE)
    #write(t(all.frame),file=filename,sep="",ncolumns = m+1,append=TRUE)
}
