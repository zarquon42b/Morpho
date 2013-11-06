#' Perform PCA based of the group means' covariance matrix
#' 
#' Calculate covariance matrix of the groupmeans and project all observations
#' into the eigenspace of this covariance matrix. This displays a low
#' dimensional between group structure of a high dimensional problem.
#' 
#' 
#' @param dataarray Either a k x m x n real array, where k is the number of
#' points, m is the number of dimensions, and n is the sample size. Or
#' alternatively a n x m Matrix where n is the numeber of observations and m
#' the number of variables (this can be PC scores for example)
#' @param groups a character/factor vector containgin grouping variable.
#' @param rounds integer: number of permutations if a permutation test of the
#' euclidean distance between group means is requested.If rounds = 0, no test
#' is performed.
#' @param tol threshold to ignore eigenvalues of the covariance matrix.
#' @param cv logical: requests leaving-one-out crossvalidation
#' @param mc.cores integer: how many cores of the Computer are allowed to be
#' used. Default is use autodetection by using detectCores() from the parallel
#' package. Parallel processing is disabled on Windows due to occasional
#' errors.
#' @param weighting logical:weight between groups covariance matrix according
#' to group sizes.
#' @return
#' \item{eigenvalues }{Non-zero eigenvalues of the groupmean covariance
#' matrix}
#' \item{groupPCs }{PC-axes - i.e. eigenvectors of the groupmean covariance
#' matrix}
#' \item{Variance }{table displaying the variance explained by eache
#' eigenvalue}
#' \item{Scores }{Scores of all observation in the PC-space}
#' \item{probs }{p-values of pairwise groupdifferences - based on
#' permuation testing}
#' \item{groupdists }{Euclidean distances between groups' averages}
#' \item{groupmeans }{Groupmeans}
#' \item{Grandmean }{Grand mean}
#' \item{CV }{Cross-validated scores}
#' @author Stefan Schlager
#' @seealso \code{\link{CVA}}
#' @references Mitteroecker P, Bookstein F 2011. Linear Discrimination,
#' Ordination, and the Visualization of Selection Gradients in Modern
#' Morphometrics. Evolutionary Biology 38:100-114.
#' 
#' Boulesteix, A. L. 2005: A note on between-group PCA, International Journal
#' of Pure and Applied Mathematics 19, 359-366.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(car)
#' data(iris)
#' vari <- iris[,1:4]
#' facto <- iris[,5]
#' pca.1 <-groupPCA(vari,groups=facto,rounds=100,mc.cores=1)
#' 
#' ### plot scores
#' scatterplotMatrix(pca.1$Scores,groups=facto, ellipse=TRUE,
#'         by.groups=TRUE,var.labels=c("PC1","PC2","PC3"))
#' 
#' ## example with shape data
#' data(boneData)
#' proc <- procSym(boneLM)
#' pop_sex <- name2factor(boneLM, which=3:4)
#' gpca <- groupPCA(proc$orpdata, groups=pop_sex, rounds=0, mc.cores=2)
#' \dontrun{
#' ## visualize shape associated with first between group PC
#' dims <- dim(proc$mshape)
#' ## calculate matrix containing landmarks of grandmean
#' grandmean <- matrix(gpca$Grandmean, dims[1], dims[2])
#' ## calculate landmarks from first between-group PC
#' #                   (+2 and -2 standard deviations)
#' gpcavis2sd<- showPC(2*sd(gpca$Scores[,1]), gpca$groupPCs, grandmean)
#' gpcavis2sd.neg<- showPC(-2*sd(gpca$Scores[,1]), gpca$groupPCs, grandmean)
#' deformGrid3d(gpcavis2sd, gpcavis2sd.neg, ngrid = 0)
#' require(rgl)
#' ## visualize grandmean mesh
#' 
#' grandm.mesh <- warp.mesh(skull_0144_ch_fe.mesh, boneLM[,,1],grandmean)
#' wire3d(grandm.mesh, col="white")
#' spheres3d(grandmean, radius=0.005)
#' }
#' 
#' 
#' @export groupPCA
groupPCA <- function(dataarray, groups, rounds = 10000,tol=1e-10,cv=TRUE,mc.cores=detectCores(), weighting=TRUE)
{
    win <- FALSE
    if(.Platform$OS.type == "windows")
        win <- TRUE
    else
        registerDoParallel(cores=mc.cores)### register parallel backend
    
    pmatrix.proc <- NULL
    proc.distout <- NULL
    lev <- NULL	
    
    groups <- factor(groups)
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
    if (1 %in% gsizes) {
        cv <- FALSE
        warning("group with one entry found - crossvalidation will be disabled.")
    }
    N <- dataarray
    if (length(dim(N)) == 3) 
        N <- vecx(N)
    N <- as.matrix(N)
    n <- dim(N)[1]
    l <- dim(N)[2]
    if (length(groups) != n)
        warning("group affinity and sample size not corresponding!")
    
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng) {
        if(gsizes[i] > 1)
            Gmeans[i, ] <- apply(N[groups==lev[i], ], 2, mean)
        else
            Gmeans[i, ] <- N[groups==lev[i], ]
    }
    if (weighting==TRUE)
        wt <- gsizes
    else
        wt <- rep(1,ng)
    wcov <- cov.wt(Gmeans,wt=wt)
    Grandm <- wcov$center
    eigenGmeans <- eigen(wcov$cov)
                                        #resGmeans <- sweep(Gmeans, 2, Grandm)
    Tmatrix <- N
    N <- sweep(N, 2, Grandm)
    valScores <- which(eigenGmeans$values > tol)
    groupScores <- N%*%(eigenGmeans$vectors[,valScores])
    groupPCs <- eigenGmeans$vectors[,valScores]

###### create a neat variance table for the groupmean PCA ###### 
    values <- eigenGmeans$values[valScores]
    if (length(values) == 1) {
        Var <- values
    } else {
        Var <- matrix(NA,length(values),3)
        Var[,1] <- values
        for (i in 1:length(values)) 
            Var[i,2] <- (values[i]/sum(values))*100
        Var[1,3] <- Var[1,2]
        for (i in 2:length(values))
            Var[i,3] <- Var[i,2]+ Var[i-1,3]
        colnames(Var) <- c("eigenvalues","% Variance","Cumulative %")
    }
### calculate between group distances ###
    proc.disto <- matrix(0, ng, ng)
    if(!is.null(lev)) {
        rownames(proc.disto) <- lev
        colnames(proc.disto) <- lev
    }	
    for (j1 in 1:(ng - 1)) 
        for (j2 in (j1 + 1):ng) 
            proc.disto[j2, j1] <- sqrt(sum((Gmeans[j1, ]- Gmeans[j2,])^2))
    
    proc.distout <- as.dist(proc.disto)

### Permutation Test for Distances	
    if (rounds > 0) {
        pmatrix.proc <- matrix(NA, ng, ng) ### generate distance matrix Euclidean
        if(!is.null(lev)) {
            rownames(pmatrix.proc) <- lev
            colnames(pmatrix.proc) <- lev
        }
        rounproc <- function(i)
            {
                shake <- sample(groups)
                dist.mat <- matrix(0,ng,ng)
                Gmeans.tmp <- matrix(0, ng, l)
                for (j in 1:ng) {
                    if(gsizes[j] > 1)
                        Gmeans.tmp[j, ] <- apply(N[shake==lev[j], ], 2, mean)
                    else
                        Gmeans.tmp[j, ] <- N[shake==lev[j], ]
                }
                dist.mat <- as.matrix(dist(Gmeans.tmp))
                return(dist.mat)
            }
        
        dist.mat.proc <- array(0, dim = c(ng, ng, rounds))
        if(win)
            a.list <- foreach(i=1:rounds)%do%rounproc(i)
        else
            a.list <- foreach(i=1:rounds)%dopar%rounproc(i)
        
        for (i in 1:rounds)
            dist.mat.proc[,,i] <- a.list[[i]]
        
        for (j1 in 1:(ng - 1)) {
            for (j2 in (j1 + 1):ng) {
                sorti <- sort(dist.mat.proc[j2, j1,])
                if (max(sorti) < proc.disto[j2, j1]) {
                    pmatrix.proc[j2, j1] <- 1/rounds
                } else {
                    marg <- min(which(sorti >= proc.disto[j2, j1]))
                    pmatrix.proc[j2, j1] <- (rounds - (marg-1))/rounds
                }
            }
        }
        pmatrix.proc <- as.dist(pmatrix.proc)
    }
    crovafun <- function(x)
        {
            crovtmp <- .groupPCAcrova(Tmatrix[-x,],groups[-x],tol=tol,groupPCs=groupPCs,weighting=weighting)
            out <- as.vector(Tmatrix[x,]-crovtmp$Grandmean) %*% as.matrix(crovtmp$PCs)
            return(out)
        }
    
    CV=NULL
    if (cv) {
        if (win)
            crossval <- foreach(i=1:n) %do% crovafun(i)
        else
            crossval <- foreach(i = 1:n) %dopar% crovafun(i)
        CV <- groupScores
        for (i in 1:n) {
            if (is.matrix(CV))
                CV[i,] <- crossval[[i]]
            else
                CV[i] <- crossval[[i]]
        }
    }
    return(list(eigenvalues=values,groupPCs=eigenGmeans$vectors[,valScores],Variance=Var,Scores=groupScores,probs=pmatrix.proc,groupdists=proc.distout,groupmeans=Gmeans,Grandmean=Grandm,CV=CV))
}
