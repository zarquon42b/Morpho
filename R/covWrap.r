covWrap <- function(s1, s2) {
    out <- .Call("covWrap", s1, s2)
    return(out)
}
covPCAwrap <- function(data,groups,rounds=1000,scramble=10) {
    out <- list()
    if (! is.factor(groups))
        groups <- as.factor(groups)
    
    groups <- factor(groups)
    lev <- levels(groups)
    nlev <- length(lev)
    for (i in 1:nlev) # center data per group
        data[groups==lev[i],] <- sweep(data[groups==lev[i],],2, apply(data[groups==lev[i],],2,mean))
    data <- as.matrix(data)
    groups <- as.integer(groups)
    out <- .Call("covPCAwrap", data, groups,scramble,rounds)
    return(out)
}
