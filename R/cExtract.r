cExtract <- function(pts.file)
{
    if (is.character(pts.file))
        x <- read.pts(pts.file)
    else
	x <- pts.file	
    
    allnames <- row.names(x)
    cs <- grep("C",allnames)
    ps <- grep("P",allnames)
    S <- grep("S",allnames)
    if (length(ps) > 0)
        cs <- c(cs,ps)
    cnames <- row.names(x)[cs]	
    olevels <- levels(as.factor(substr(cnames,1,4)))
    S <- grep("S",allnames)
    if (length(S) == 0 ) {
        S <- NULL
    } else {
        S <- "S"
    }
    olevels <- c(S, olevels)
    tl <- length(olevels)
    
    out <- list()	
    for (i in 1:tl)
        out[[olevels[i]]] <- grep(olevels[i],allnames)
    
    return(out)
}
