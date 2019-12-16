#' @rdname mesh2ply
#' @importFrom Rvcg vcgObjWrite
#' @export
mesh2obj <- function(x,filename=dataname, writeNormals = TRUE ) {
    dataname <- deparse(substitute(x))
    vcgObjWrite(x,filename=dataname,writeNormals)
    
}
## {
##     ismatrix <- FALSE
##     x.it <- NULL
##     dataname <- deparse(substitute(x))
##     if (is.matrix(x)) {
##         ismatrix <- TRUE
##         dimsx <- dim(x)
##         if (dimsx[2] == 3 && dimsx[1] != 3)
##             x <- t(x)
##         x <- list(vb=x)
##     }        
##     x.vb <- cbind("v",t(x$vb[1:3,]))
##     if (! ismatrix)
##         x.it <- cbind("f",t(x$it))
##     obj <- rbind(x.vb,x.it)
    
##     write.obj(obj, filename=filename)
## }

