#' calculate typicality probabilities
#' 
#' get the probability for an observation belonging to a given multivariate
#' nromal distribution
#' 
#' 
#' @title calculate typicality probabilities
#' @param x vector or matrix of data the probability is to be calculated.
#' @param data Reference data set.
#' @param groups vector containing grouping information.
#' @param small adjustion of Mahalanobis D^2 for small sample sizes as
#' suggested by Wilson (1981), only takes effect when method="wilson".
#' @param method select method for probability estimation. Available options
#' are "chisquare" (or any abbreviation) or "wilson". "chisquare" exploits
#' simply the chisquare distribution of the mahalanobisdistance, while "wilson"
#' uses the methods suggested by Wilson(1981). Results will be similar in
#' general.
#' @param outlier probability threshold below which a specimen will not be
#' assigned to any group-
#' @param sep logical: if TRUE, probability will be calculated from the pooled
#' within group covariance matrix.
#' @param center vector: specify custom vector to calculate distance to. If not
#' defined, group mean will be used.
#' @param cova covariance matrix to calculate mahalanobis-distance: specify
#' custom covariance matrix to calculate distance.
#' @return typprob: returns a vector of probabilities.
#' 
#' typprobClass:
#' \item{probs }{matrix of probabilities for each group}
#' \item{groupaffin }{vector of groups each specimen has been assigned to.
#' Outliers are classified "none"}
#' @author Stefan Schlager
#' @references Albrecht G. 1992. Assessing the affinities of fossils using
#' canonical variates and generalized distances Human Evolution 7:49-69.
#' 
#' Wilson S. 1981. On comparing fossil specimens with population samples
#' Journal of Human Evolution 10:207 - 214.
#' 
#' @examples
#' 
#' library(shapes)
#' data <- procSym(gorf.dat)$PCscores[,1:3]
#' probas <- typprob(data,data,small=TRUE)### get probability for each specimen
#' 
#' ### now we check how this behaves compared to the mahalanobis distance
#' maha <- mahalanobis(data,apply(data,2,mean),cov(data))
#' plot(probas,maha,xlab="Probability",ylab="Mahalanobis D^2")
#' 
#' data2 <- procSym(abind(gorf.dat,gorm.dat))$PCscores[,1:3]
#' fac <- as.factor(c(rep("female",dim(gorf.dat)[3]),rep("male",dim(gorm.dat)[3])))
#' typClass <- typprobClass(data2,data2,fac,method="w",small=TRUE)
#' ## only 59 specimen is rather small.
#' typClass2 <- typprobClass(data2,data2,fac,method="c")## use default settings
#' 
#' ### check results for first method:
#' ct <- table(fac,typClass$groupaffin)
#' ct #view classification table
#' ### get percentage of correct classification
#' prop.table(ct, 1)
#' 
#' ### check results for second method:
#' ct1 <- table(fac,typClass2$groupaffin)
#' ct1 #view classification table ### one specimen has been tagged an outlier.
#' ### get percentage of correct callification
#' prop.table(ct1, 1) 
#' 
#' 
#' @rdname typprob
#' @export
typprob <- function(x,data,small=FALSE, method=c("chisquare","wilson"), center=NULL, cova=NULL)
{
    method <- substr(method[1],1L,1L)
    if (is.matrix(x))
        nx <- dim(x)[2]
    else
        nx <- length(x)

    ndata <- dim(data)[1]
    if (is.null(center))
        center <- apply(data,2,mean)
    
    if (is.null(cova))
        cova <- cov(data)
    
    dists <- mahalanobis(x, center=center, cov=cova)
    if (method == "w")
        {
            if (small)
                dists <- ndata*log(1+(dists*ndata/(ndata^2-1)))
            t2 <- ndata*dists/(ndata+1)
            f <- t2*(ndata-nx)/(nx*(ndata-1))
            alpha <- pf(f,nx,(ndata-nx),lower.tail=T)
            alpha <- 1-alpha
        }
    else if (method == "c")
        alpha <- pchisq(dists,nx,lower.tail=F)
    else
        stop("please chose valid method")
    
    return(alpha)
}

#' @rdname typprob
#' @export
typprobClass <- function(x,data,groups,small=FALSE,method=c("chisquare","wilson"),outlier=0.01,sep=FALSE)
{
    if (!is.factor(groups))
        {
            groups <- as.factor(groups)
            warning("groups coerced to factors")
        }
    probs <- NULL
    groups <- factor(groups)
    glev <- levels(groups)
    nlev <- length(glev)
    cova <- NULL
    if (sep)
        cova <- covW(data,groups)
    for( i in 1:nlev)
        {
            if (sep)# use separate covariance matrices for scaling
                tmp <- typprob(x,data[groups==glev[i],],small=small,method=method)
            else #used pooled within group covariance
                tmp <- typprob(x,data,small=small,method=method,center=apply(data[groups==glev[i],],2,mean),cova=cova)
            probs <- cbind(probs,tmp)
        }
    colnames(probs) <- as.character(glev)
    classify <- apply(probs,1,function(x){out <- which(x==max(x));return(out)})
    glev <- c(glev,"none")
    outsider <- which(apply(probs,1,max) < outlier)
    classify[outsider] <- nlev+1
    groupaffin <- as.factor(glev[classify])
    
    out <- (list(probs=probs,groupaffin=groupaffin))
    class(out) <- "typprob"
    return(out)
}

#' calculate the pooled within groups covariance matrix
#'
#' calculate the pooled within groups covariance matrix
#' @param data a matrix containing data
#' @param groups grouping variables
#' @return Returns the pooled within group covariance matrix
#' @author Stefan Schlager
#' @seealso \code{\link{cov}}, \code{\link{typprobClass}}
#' @examples
#' data(iris)
#' poolCov <- covW(iris[,1:4],iris[,5])
#' @export
covW <- function(data, groups)
    {
        if (!is.factor(groups)) {
            groups <- as.factor(groups)
            warning("groups coerced to factors")
        }
        ndata <- dim(data)[1]
        groups <- factor(groups)
        glev <- levels(groups)
        nlev <- length(glev)
        gsizes <- as.vector(tapply(groups, groups, length))
        covWithin <- 0
        for (i in 1:nlev)
            if (gsizes[i] > 1) 
                covWithin <- covWithin + (cov(data[groups==glev[i],]) * (gsizes[i]-1))
        
        covWithin <- covWithin/(ndata - nlev)
        return(covWithin)
    }
