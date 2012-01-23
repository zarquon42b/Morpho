addo<-function(a3) 
{	## copy of add in shapes package
    s <- 0
    for (i in 1:dim(a3)[3]) {
        s <- s + a3[, , i]
    }
    return(s)
}

