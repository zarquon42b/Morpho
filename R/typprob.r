#' calculate typicality probabilities
#' 
#' get the probability for an observation belonging to a given multivariate
#' nromal distribution
#' 
#' 
#' @title calculate typicality probabilities
#' @param x vector or matrix of data the probability is to be calculated.
#' @param data Reference data set. If missing x will be used.
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
#' @param cv logical: if data is missing and \code{cv=TRUE}, the resulting classification will be validated by leaving-one-out crossvalidation.
#' @return typprob: returns a vector of probabilities.
#' 
#' typprobClass:
#' \item{probs }{matrix of probabilities for each group}
#' \item{groupaffin }{vector of groups each specimen has been assigned to. Outliers are classified "none"}
#' \item{probsCV }{cross-validated matrix of probabilities for each group}
#' \item{groupaffinCV }{cross-validated vector of groups each specimen has been assigned to. Outliers are classified "none"}
#' \item{self }{logical: if TRUE, the data has been classified by self-inference.}
#'
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
#' maha <- mahalanobis(data,colMeans(data),cov(data))
#' plot(probas,maha,xlab="Probability",ylab="Mahalanobis D^2")
#' 
#' data2 <- procSym(abind(gorf.dat,gorm.dat))$PCscores[,1:3]
#' fac <- as.factor(c(rep("female",dim(gorf.dat)[3]),rep("male",dim(gorm.dat)[3])))
#' typClass <- typprobClass(data2,groups=fac,method="w",small=TRUE,cv=TRUE)
#' ## only 59 specimen is rather small.
#' typClass2 <- typprobClass(data2,groups=fac,method="c",cv=TRUE)## use default settings
#' 
#' ### check results for first method:
#' typClass
#' 
#' 
#' ### check results for second method:
#' typClass2
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
        center <- colMeans(data)
    
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
typprobClass <- function(x,data,groups,small=FALSE,method=c("chisquare","wilson"),outlier=0.01,sep=FALSE,cv=TRUE)
{
    if (!is.factor(groups)) {
            groups <- as.factor(groups)
            warning("groups coerced to factors")
        }
    self <- FALSE
    if (missing(data)) {
        data <- x
        self <- TRUE
    }
    probs <- NULL
    groups <- factor(groups)
    glev <- levels(groups)
    nlev <- length(glev)
    cova <- NULL
    if (!sep)
        cova <- covW(data,groups)
    for( i in 1:nlev)
        {
            if (sep)# use separate covariance matrices for scaling
                tmp <- typprob(x,data[groups==glev[i],,drop=FALSE],small=small,method=method)
            else #used pooled within group covariance
                tmp <- typprob(x,data,small=small,method=method,center=colMeans(data[groups==glev[i],,drop=FALSE]),cova=cova)
            probs <- cbind(probs,tmp)
        }
    colnames(probs) <- as.character(glev)
    classify <- apply(probs,1,function(x){out <- which(x==max(x));return(out)})
    glev <- c(glev,"none")
    outsider <- which(apply(probs,1,max) < outlier)
    classify[outsider] <- nlev+1
    groupaffin <- as.factor(glev[classify])
    CV <- probsCV <- gaffinCV <- NULL
    if (cv && self) {
        probsCV <- probs
        gaffinCV <- as.character(groupaffin)
        
        for(i in 1:length(groups)) {
            tmp <- typprobClass(data[i,],data[-i,,drop=FALSE],groups=groups[-i],small=small,method=method,outlier = outlier, sep=sep,cv=FALSE)
            probsCV[i,] <- tmp$probs
            gaffinCV[i] <- as.character(tmp$groupaffin)
        }
    } else
        cv <- FALSE
    out <- (list(probs=probs,groupaffin=groupaffin,cv=cv,self=self))
    out$probsCV <- probsCV
    out$groupaffinCV <- factor(gaffinCV)
    out$groups <- groups
    
    class(out) <- "typprob"
    return(out)
}

#' @export
print.typprob <- function(x,...) {
    print(classify(x,cv=TRUE))
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
                covWithin <- covWithin + (cov(data[groups==glev[i],,drop=FALSE]) * (gsizes[i]-1))
        
        covWithin <- covWithin/(ndata - nlev)
        return(covWithin)
    }
