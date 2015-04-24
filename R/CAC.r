#' calculate common allometric component
#'
#' calculate common allometric component
#'
#' @param x datamatrix (e.g. with PC-scores) or 3D-array with landmark coordinates
#' @param size vector with Centroid sizes
#' @param groups grouping variable
#' @param log logical: use \code{log(size)}
#'
#' @return
#' \item{CACscores}{common allometric component scores}
#' \item{CAC}{common allometric component}
#' \item{x}{(group-) centered data}
#' \item{sc}{CAC reprojected into original space by applying \code{CAC \%*\% x}}
#' \item{RSCscores}{residual shape component scores}
#' \item{RSC}{residual shape components}
#' \item{gmeans}{groupmeans}
#' \item{CS}{the centroid sizes (log transformed if \code{log = TRUE})}
#'
#' @references
#' Mitteroecker P, Gunz P, Bernhard M, Schaefer K, Bookstein FL. 2004. Comparison of cranial ontogenetic trajectories among great apes and humans. Journal of Human Evolution 46(6):679-97.
#' @examples
#' data(boneData)
#' proc <- procSym(boneLM)
#' pop.sex <- name2factor(boneLM,which=3:4)
#' cac <- CAC(proc$rotated,proc$size,pop.sex)
#' plot(cac$CACscores,cac$size)#plot scores against Centroid size
#' cor.test(cac$CACscores,cac$size)#check for correlation
#' #visualize differences between large and small on the sample's consensus
#' \dontrun{
#' large <- showPC(max(cac$CACscores),cac$CAC,proc$mshape)
#' small <- showPC(min(cac$CACscores),cac$CAC,proc$mshape)
#' deformGrid3d(small,large,ngrid=0)
#' }
#' @export

CAC <- function(x, size, groups=NULL, log=FALSE) {
    if (length(dim(x)) == 3)
        x <- vecx(x)
    if (!is.null(groups)) {
        groups <- factor(groups)
        lm0 <- lm(x ~ groups)
        x <- lm0$residuals
        lev <- levels(groups)
        nlev <- length(lev)
        indices <- sapply(1:nlev,function(x) x <- which(groups==lev[x])[1])
        gmeans <- lm0$fitted.values[indices,]
        rownames(gmeans) <- lev
    } else {
        x <- scale(x,scale = FALSE)
        gmeans <- colMeans(x)
    }
    if (log)
        size <- log(size)
    allo <- (t(x)%*%size)/sqrt(sum(size^2))
    allo <- allo/(sqrt(sum(allo^2)))
    xout <- x%*%allo
    shapechange <- xout%*%t(allo)
    W <- x-shapechange
    WTW <- svd(crossprod(W))
    RSC <- WTW$u[,which(WTW$d > 1e-13)]
    RSCscores <- x%*%RSC
    
    return(list(CACscores=xout,CAC=allo,x=x,sc=shapechange,RSCscores=RSCscores, RSC=RSC,gmeans=gmeans,size=size))
}
