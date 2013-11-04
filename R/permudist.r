permudist <- function(data, groups, rounds=1000, which=1:2, mc.cores = detectCores())
{
    win <- FALSE
    if(.Platform$OS.type == "windows")
        win <- TRUE
    else
        registerDoParallel(cores=mc.cores)### register parallel backend
    out <- list()
### configure grouping ####
    N <- data
    if (is.vector(N)) {
        N <- as.matrix(N)
    } else if (dim(N)[2] == 3)
        N <- vecx(N)
    
    if (!is.factor(groups))
        groups <- factor(groups)
    
    if (is.factor(groups)) {
        groups <- factor(groups)
        lev <- levels(groups)
        #levn <- length(lev)
    }
    old_groups <- groups
    groups <- factor(groups[groups %in% lev[which]])
    lev <- levels(groups)
    N <- N[which(old_groups %in% lev),]

### end configure grouping ####
    #n <- dim(N)[1]
   
    if (dim(N)[2] == 1) {
        mean1 <- mean(N[groups == lev[1]])
        mean2 <- mean(N[groups == lev[2]])
    } else {
        mean1 <- apply(N[groups == lev[1],], 2, mean)
        mean2 <- apply(N[groups == lev[2],], 2, mean)
    }
    dist <- sqrt(sum((mean1-mean2)^2))
    
    if (rounds > 0) {
        permu <- function(x)
            {
                shake <- sample(groups)
                disto <- permudist(N, shake, rounds = 0, mc.cores=1)$dist
                return(disto)
            }
        i <- 0
        if (win)
            dists <- foreach(i = 1:rounds,.combine=c) %do% permu(i)
        else
            dists <- foreach(i = 1:rounds,.combine=c) %dopar% permu(i)
        
        p.value <- length(which(dists >= dist))
        if (p.value > 0) {
            p.value <- p.value/rounds
            names(p.value) <- "p-value"
        } else {
            p.value <- 1/rounds
            names(p.value) <- "p-value <"
        }
        out$p.value <- p.value
        out$permudist <- dists
    }
    
    out$dist <- dist

return(out)
}
