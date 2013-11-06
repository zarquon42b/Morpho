#' performs permutation testing for group differences.
#' 
#' This function compares the distance between two groupmeans to the distances
#' obtained by random assignment of observations to this groups.
#' 
#' 
#' @param data array or matrix containing data
#' @param groups factors determining grouping.
#' @param rounds number of permutations
#' @param which in case the factor levels are >2 this determins which
#' factorlevels to use
#' @param mc.cores integer: determines how many cores to use for the
#' computation. The default is autodetect. But in case, it doesn't work as
#' expected cores can be set manually. Parallel processing is disabled on
#' Windows due to occasional errors.
#' @return
#' \item{permudist }{vector containing distances between random group means}
#' \item{dist }{distance between actual group means}
#' \item{p.value }{p-value obtained by comparing the actual distance to randomly acquired distances}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' perm <- permudist(proc$PCscores[,1:10], groups=groups, mc.cores=2, rounds=100)
#' \dontrun{
#' #visualize results
#' hist(perm$permudist, xlim=c(0,0.1),main="measured vs. random distances",
#'      xlab="distances")#random distances
#' points(perm$dist,50,col=2,pch=19)#actual distance
#' text(perm$dist,70,label=paste("actual distance\n(p=",perm$p.value,")"))
#' 
#' ## now we concentrate only on sex dimorphism between Europeans
#' groups <- name2factor(boneLM,which=3:4)
#' levels(groups)
#' perm1 <- permudist(proc$PCscores, groups=groups,which=3:4, mc.cores=2, rounds=100)
#' # plot histogram of random distances
#' hist(perm1$permudist, xlim=c(0,0.1),main="measured vs. random distances",
#'      xlab="distances")
#' points(perm1$dist,10,col=2,pch=19)#actual distance
#' text(perm1$dist,15,label=paste("actual
#' distance\n(p=",perm1$p.value,")"))
#' }
#' 
#' @export permudist
permudist <- function(data, groups, rounds=1000, which=1:2, mc.cores = detectCores())
{
    win <- FALSE
    if(.Platform$OS.type == "windows" || mc.cores == 1)
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
