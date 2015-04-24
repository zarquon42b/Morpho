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
#' \item{Variance }{table displaying the between-group variance explained by each between group PC}
#' \item{Scores }{Scores of all observation in the PC-space}
#' \item{probs }{p-values of pairwise groupdifferences - based on
#' permuation testing}
#' \item{groupdists }{Euclidean distances between groups' averages}
#' \item{groupmeans }{Groupmeans}
#' \item{Grandmean }{Grand mean}
#' \item{CV }{Cross-validated scores}
#' \item{groups }{grouping Variable}
#' \item{resPCs}{PCs orthogonal to the between-group PCs}
#' \item{resPCscores}{Scores of the residualPCs}
#' \item{resVar}{table displaying the residual variance explained by each residual PC}
#' \item{combinedVar}{table displaying the overall variance explained by the between-group PCs and residual PC. Check the rownames to identify which type belongs to which value}
#' @author Stefan Schlager
#' @seealso \code{\link{CVA}}
#' @references Mitteroecker P, Bookstein F 2011. Linear Discrimination,
#' Ordination, and the Visualization of Selection Gradients in Modern
#' Morphometrics. Evolutionary Biology 38:100-114.
#' 
#' Boulesteix, A. L. 2005: A note on between-group PCA, International Journal
#' of Pure and Applied Mathematics 19, 359-366.
#' 
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
#' gpcavis2sd<- showPC(2*sd(gpca$Scores[,1]), gpca$groupPCs[,1], grandmean)
#' gpcavis2sd.neg<- showPC(-2*sd(gpca$Scores[,1]), gpca$groupPCs[,1], grandmean)
#' deformGrid3d(gpcavis2sd, gpcavis2sd.neg, ngrid = 0)
#' require(rgl)
#' ## visualize grandmean mesh
#' 
#' grandm.mesh <- tps3d(skull_0144_ch_fe.mesh, boneLM[,,1],grandmean)
#' wire3d(grandm.mesh, col="white")
#' spheres3d(grandmean, radius=0.005)
#' }
#' 
#' 
#' @export
groupPCA <- function(dataarray, groups, rounds = 10000,tol=1e-10,cv=TRUE,mc.cores=parallel::detectCores(), weighting=TRUE)
{
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
    rownames(Gmeans) <- lev
    for (i in 1:ng) {
        Gmeans[i, ] <- colMeans(N[groups==lev[i], ,drop=F])
    }
    if (weighting)
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
    residuals <- N-groupScores%*%t(groupPCs)
    resPrcomp <- prcomp(residuals,center = F,tol=sqrt(tol))
   
    
###### create a neat variance table for the groupmean PCA ###### 
    values <- eigenGmeans$values[valScores]
    bgnames <- c(paste("bgPC",1:length(values),sep="_"))
    Var <- createVarTable(values,FALSE,rownames = bgnames)
    cnames <- c(paste("bgPC",1:length(values),sep="_"),paste("resPC",1:length(resPrcomp$sdev),sep="_"))
    combinedVar <- createVarTable(c(values,resPrcomp$sdev^2),square = FALSE,rownames = cnames)
    resnames <- paste("resPC",1:length(resPrcomp$sdev),sep="_")
    resVar <- createVarTable(resPrcomp$sdev,square = TRUE,rownames = resnames)
### calculate between group distances ###
    
### Permutation Test for Distances	
    
    shakeitall <- permudist(N,groups,rounds=rounds)
    if (rounds > 0)
        pmatrix.proc <- shakeitall$p.value
    else
        pmatrix.proc <- NULL
    proc.distout <- shakeitall$dist
    
    crovafun <- function(x)
        {
        crovtmp <- .groupPCAcrova(Tmatrix[-x,,drop=F],groups[-x],tol=tol,groupPCs=groupPCs,weighting=weighting)
        out <- as.vector(Tmatrix[x,]-crovtmp$Grandmean) %*% as.matrix(crovtmp$PCs)
        return(out)
    }

CV=NULL
    if (cv) {
        win <- FALSE
        if(.Platform$OS.type == "windows")
            win <- TRUE
        else
            registerDoParallel(cores=mc.cores)### register parallel backend
        
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
    out <- list(eigenvalues=values,groupPCs=eigenGmeans$vectors[,valScores],Variance=Var,Scores=groupScores,probs=pmatrix.proc,groupdists=proc.distout,groupmeans=Gmeans,Grandmean=Grandm,CV=CV,groups=groups,resPCs=resPrcomp$rotation,resPCscores=resPrcomp$x,resVar=resVar,combinedVar=combinedVar)
    class(out) <- "bgPCA"
    return(out)
}


