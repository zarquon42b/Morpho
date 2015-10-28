#' Canonical Variate Analysis
#' 
#' performs a Canonical Variate Analysis.
#' 
#' 
#' @param dataarray Either a k x m x n real array, where k is the number of
#' points, m is the number of dimensions, and n is the sample size. Or
#' alternatively a n x m Matrix where n is the numeber of observations and m
#' the number of variables (this can be PC scores for example)
#' @param groups a character/factor vector containgin grouping variable.
#' @param weighting Logical: Determines whether the between group covariance 
#' matrix and Grandmean is to be weighted according to group size.
#' @param tolinv Threshold for the eigenvalues of the pooled
#' within-group-covariance matrix to be taken as zero - for calculating the
#' general inverse of the pooled withing groups covariance matrix.
#' @param plot Logical: determins whether in the two-sample case a histogramm
#' ist to be plotted.
#' @param rounds integer: number of permutations if a permutation test of the
#' Mahalanobis distances (from the pooled within-group covariance matrix) and Euclidean distance between group means is requested.If
#' rounds = 0, no test is performed.
#' @param cv logical: requests a Jackknife Crossvalidation.
#' @param p.adjust.method method to adjust p-values for multiple comparisons see \code{\link{p.adjust.methods}} for options.
#' @return
#' \item{CV }{A matrix containing the Canonical Variates}
#' \item{CVscores }{A matrix containing the individual Canonical Variate scores}
#' \item{Grandm }{a vector or a matrix containing the Grand Mean (depending if the input is an array or a matrix)}
#' \item{groupmeans }{a matrix or an array containing the group means (depending if the input is an array or a matrix)}
#' \item{Var }{Variance explained by the Canonical Variates}
#' \item{CVvis }{Canonical Variates projected back into the original space - to be used for visualization purposes, for details see example below}
#' \item{Dist }{Mahalanobis Distances between group means - if requested
#' tested by permutation test if the input is an array it is assumed to be
#' superimposed Landmark Data and Procrustes Distance will be calculated}
#' \item{CVcv }{A matrix containing crossvalidated CV scores}
#' \item{groups }{factor containing the grouping variable}
#' \item{class }{classification results based on posteriror probabilities. If cv=TRUE, this will be done by a leaving-one-out procedure}
#'  \item{posterior }{posterior probabilities}
#'  \item{prior }{prior probabilities}
#' @author Stefan Schlager
#' @seealso \code{\link{groupPCA}}
#' @references Cambell, N. A. & Atchley, W. R.. 1981 The Geometry of Canonical
#' Variate Analysis: Syst. Zool., 30(3), 268-280.
#' 
#' Klingenberg, C. P. & Monteiro, L. R. 2005 Distances and directions in
#' multidimensional shape spaces: implications for morphometric applications.
#' Systematic Biology 54, 678-688.
#' @examples
#' 
#' ## all examples are kindly provided by Marta Rufino
#' 
#' library(shapes)
#' # perform procrustes fit on raw data
#' alldat<-procSym(abind(gorf.dat,gorm.dat))
#' # create factors
#' groups<-as.factor(c(rep("female",30),rep("male",29)))
#' # perform CVA and test Mahalanobis distance
#' # between groups with permutation test by 100 rounds)            
#' cvall<-CVA(alldat$orpdata,groups,rounds=10000)     
#' ## visualize a shape change from score -5 to 5:
#' cvvis5 <- 5*matrix(cvall$CVvis[,1],nrow(cvall$Grandm),ncol(cvall$Grandm))+cvall$Grandm
#' cvvisNeg5 <- -5*matrix(cvall$CVvis[,1],nrow(cvall$Grandm),ncol(cvall$Grandm))+cvall$Grandm
#' plot(cvvis5,asp=1)
#' points(cvvisNeg5,col=2)
#' for (i in 1:nrow(cvvisNeg5))
#'   lines(rbind(cvvis5[i,],cvvisNeg5[i,]))
#' 
#' ### Morpho CVA
#' data(iris)
#' vari <- iris[,1:4]
#' facto <- iris[,5]
#' 
#' cva.1=CVA(vari, groups=facto)
#' ## get the typicality probabilities and resulting classifications - tagging
#' ## all specimens with a probability of < 0.01 as outliers (assigned to no class)
#' typprobs <- typprobClass(cva.1$CVscores,groups=facto)
#' print(typprobs)
#' ## visualize the CV scores by their groups estimated from (cross-validated)
#' ## typicality probabilities:
#' require(car)
#' scatterplot(cva.1$CVscores[,1],cva.1$CVscores[,2],groups=typprobs$groupaffinCV,
#'                   smooth=FALSE,reg.line=FALSE)
#' # plot the CVA
#' plot(cva.1$CVscores, col=facto, pch=as.numeric(facto), typ="n",asp=1,
#'    xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"%")),
#'    ylab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"%")))
#'   
#'   text(cva.1$CVscores, as.character(facto), col=as.numeric(facto), cex=.7)
#' 
#'   # add chull (merge groups)
#'   for(jj in 1:length(levels(facto))){
#'         ii=levels(facto)[jj]
#'     kk=chull(cva.1$CVscores[facto==ii,1:2])
#'     lines(cva.1$CVscores[facto==ii,1][c(kk, kk[1])],
#'     cva.1$CVscores[facto==ii,2][c(kk, kk[1])], col=jj)
#'     }
#' 
#'   # add 80% ellipses
#'   require(car)
#'   for(ii in 1:length(levels(facto))){
#'     dataEllipse(cva.1$CVscores[facto==levels(facto)[ii],1],
#'     cva.1$CVscores[facto==levels(facto)[ii],2], 
#'                     add=TRUE,levels=.80, col=c(1:7)[ii])}
#' 
#'   # histogram per group
#'   require(lattice)
#'   histogram(~cva.1$CVscores[,1]|facto,
#'   layout=c(1,length(levels(facto))),
#'           xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"%")))
#'   histogram(~cva.1$CVscores[,2]|facto, layout=c(1,length(levels(facto))),
#'           xlab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"%")))
#' 
#'   # plot Mahalahobis
#'   dendroS=hclust(cva.1$Dist$GroupdistMaha)
#'   dendroS$labels=levels(facto)
#'   par(mar=c(4,4.5,1,1))
#'   dendroS=as.dendrogram(dendroS)
#'   plot(dendroS, main='',sub='', xlab="Geographic areas",
#'           ylab='Mahalahobis distance')
#' 
#'  
#'    # Variance explained by the canonical roots:
#'    cva.1$Var
#'    # or plot it:
#'    barplot(cva.1$Var[,2])
#' 
#' # another landmark based example in 3D: 
#' data(boneData)
#' groups <- name2factor(boneLM,which=3:4)
#' proc <- procSym(boneLM)
#' cvall<-CVA(proc$orpdata,groups)    
#' #' ## visualize a shape change from score -5 to 5:
#' cvvis5 <- 5*matrix(cvall$CVvis[,1],nrow(cvall$Grandm),ncol(cvall$Grandm))+cvall$Grandm
#' cvvisNeg5 <- -5*matrix(cvall$CVvis[,1],nrow(cvall$Grandm),ncol(cvall$Grandm))+cvall$Grandm
#' \dontrun{
#' #visualize it
#' deformGrid3d(cvvis5,cvvisNeg5,ngrid = 0)
#' }
#' 
#' #for using (e.g. the first 5) PCscores, one will do:
#' cvall <- CVA(proc$PCscores[,1:5],groups)    
#' #' ## visualize a shape change from score -5 to 5:
#' cvvis5 <- 5*cvall$CVvis[,1]+cvall$Grandm
#' cvvisNeg5 <- -5*cvall$CVvis[,1]+cvall$Grandm
#' cvvis5 <- showPC(cvvis5,proc$PCs[,1:5],proc$mshape)
#' cvvisNeg5 <- showPC(cvvisNeg5,proc$PCs[,1:5],proc$mshape)
#' \dontrun{
#' #visualize it
#' deformGrid3d(cvvis5,cvvisNeg5,ngrid = 0)
#' }
#' @export
CVA <- function (dataarray, groups, weighting = TRUE, tolinv = 1e-10,plot = TRUE, rounds = 0, cv = FALSE,p.adjust.method= "none") 
{
    groups <- factor(groups)
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
    
    if (1 %in% gsizes) {
        cv <- FALSE
        warning("group with one entry found - crossvalidation will be disabled.")
    }
    prior <- gsizes/sum(gsizes)
    N <- dataarray
    n3 <- FALSE
    if (length(dim(N)) == 3) {
        k <- dim(N)[1]
        m <- dim(N)[2]
        N <- vecx(N)
        n3 <- TRUE
    }
    N <- as.matrix(N)
    n <- dim(N)[1]
    l <- dim(N)[2]
    if (length(groups) != n)
        warning("group affinity and sample size not corresponding!")
    
   
    if (is.vector(N) || dim(N)[2] == 1)
        stop("data should contain at least 2 variable dimensions")
    
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng) {
        Gmeans[i, ] <- colMeans(N[groups==lev[i], ,drop=FALSE])
    }
    if (weighting) {
        Grandm <- colSums(Gmeans*gsizes)/n ## calculate weighted Grandmean (thanks to Anne-Beatrice Dufour for the bug-fix)
    } else {
        Grandm <- colMeans(Gmeans)
    }
    Tmatrix <- N
    N <- sweep(N, 2, Grandm) #center data according to Grandmean
    resGmeans <- sweep(Gmeans, 2, Grandm)
    
    if (weighting) {
        for (i in 1:ng) 
            resGmeans[i, ] <- sqrt(gsizes[i]) * resGmeans[i, ]
        X <- resGmeans
    } else 
        X <- sqrt(n/ng) * resGmeans

    
    
    covW <- covW(N, groups)
    eigW <- eigen(covW*(n - ng))
    eigcoW <- eigen(covW)
    U <- eigW$vectors
    E <- eigW$values
    Ec <- eigcoW$values
    Ec2 <- Ec
    
    if (min(E) < tolinv)
        cat(paste("singular Covariance matrix: General inverse is used. Threshold for zero eigenvalue is", tolinv, "\n"))
    for (i in 1:length(eigW$values)) {
        if (Ec[i] < tolinv) {
            E[i] <- Ec[i] <- Ec2[i] <- 0
        } else {
            E[i] <- sqrt(1/E[i])
            Ec[i] <- sqrt(1/Ec[i])
            Ec2[i] <- (1/Ec2[i])
        }
    }
    invcW <- diag(Ec)
    irE <- diag(E)
    ZtZ <- irE %*% t(U) %*% t(X) %*% X %*% U %*% irE
    eigZ <- eigen(ZtZ,symmetric=TRUE)
    useEig <- min((ng-1), l)
    A <- Re(eigZ$vectors[, 1:useEig])
    CV <- U %*% invcW %*% A
    CVvis <- covW %*% CV
    CVscores <- N %*% CV

    roots <- eigZ$values[1:useEig]
    if (length(roots) == 1) {
        Var <- matrix(roots, 1, 1)
        colnames(Var) <- "Canonical root"
    } else {
        Var <- matrix(NA, length(roots), 3)
        Var[, 1] <- as.vector(roots)
        for (i in 1:length(roots)) 
            Var[i, 2] <- (roots[i]/sum(roots)) * 100
        Var[1, 3] <- Var[1, 2]
        for (i in 2:length(roots))
            Var[i, 3] <- Var[i, 2] + Var[i - 1, 3]
        colnames(Var) <- c("Canonical roots", "% Variance", "Cumulative %")
    }
    if (plot == TRUE && ng == 2) {
        histGroup(CVscores,groups = groups)
    }
    U2 <- eigcoW$vectors
    winv <- U2 %*% (diag(Ec2)) %*% t(U2)
    disto <- matrix(0, ng, ng)
    rownames(disto) <- colnames(disto) <- lev
         
    #proc.distout <- NULL
    #pmatrix <- NULL
    #pmatrix.proc <- NULL


### Permutation Test for Distances	
    Dist <- .CVAdists(N, groups, rounds,  winv ,p.adjust.method)

    if (n3) {
        Grandm <- matrix(Grandm, k,m)
        groupmeans <- array(as.vector(t(Gmeans)), dim = c(k,m,ng))
        dimnames(groupmeans)[[3]] <- lev
    } else {
        groupmeans <- Gmeans
        rownames(groupmeans) <- lev
    }
    classVec <- groups
    classprobs <- NULL
    classdist <- NULL
    CVcv <- NULL
    if (cv == TRUE) {
        CVcv <- CVscores
        ng <- length(groups)
        a.list <- as.list(1:n)
        crovafun <- function(i)
            {
                tmp <- .CVAcrova(Tmatrix[-i, ],groups=groups[-i],test=CV, tolinv = tolinv, weighting=weighting)
                out <- (Tmatrix[i, ]-tmp$Grandmean) %*% tmp$CV
                tmpdist <- rowSums(sweep(tmp$meanscores,2,as.vector(out))^2)
                post <- probpost(tmpdist,prior)
                myclass <- names(tmpdist)[which(post == max(post))]
                return(list(scores=out,class=myclass,dist=tmpdist,post=post))
            }
        a.list <- lapply(a.list, crovafun)
        
        for (i in 1:n) {
            CVcv[i,] <- a.list[[i]]$scores
            classVec[i] <- a.list[[i]]$class
            classprobs <- rbind(classprobs,a.list[[i]]$post)
        }
        rownames(classprobs) <- rownames(N)
        colnames(classprobs) <- lev
        names(classVec) <- rownames(N)
        #classVec <- factor(classVec)
    } 
       
    out <- list(CV = CV, CVscores = CVscores, Grandm = Grandm,
                groupmeans = groupmeans, Var = Var, CVvis = CVvis,
                Dist = Dist, CVcv = CVcv,groups=groups,class=classVec,posterior=classprobs,prior=prior
                )
    class(out) <- "CVA"
    if (!cv) {
        cl <- classify(out,cv=FALSE)
        out$class <- cl$class
        out$posterior <- cl$posterior
    }

    return(out)
}


probpost <- function(dist,prior) {
    posts <- NULL
    const <- sqrt(2*pi)
    for(i in 1:length(dist))        
        posts[i] <- const*exp(-0.5*dist[i])*prior[i]
    posts <- posts/sum(posts)
    return(posts)
}
    
#' @export
print.CVA <- function(x,...) {
    print(classify(x,cv=TRUE))
}

#' @export
print.bgPCA <- function(x,...) {
    print(classify(x,cv=TRUE))
}
