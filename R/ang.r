ang <- function(x,y,circle=TRUE) {
    a <- .Call("ang_calcC",x,y,circle)
    return(a)
}
