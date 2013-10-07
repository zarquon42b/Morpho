cExtract <- function(pts.file)
{
    if (is.character(pts.file))
        x <- read.pts(pts.file)
    else
	x <- pts.file	
    
    allnames <- row.names(x)
    cs <- grep("C",allnames)
    ps <- grep("P",allnames)
    if (length(ps)>0)
        cs <- c(cs,ps)
    cnames <- row.names(x)[cs]	
    t <- levels(as.factor(substr(cnames,1,4)))
    t <- c("S",t)
    tl <- length(t)
    
    out <- list()	
    for (i in 1:tl)
        out[[t[i]]] <- grep(t[i],allnames)
    
    return(out)
}
