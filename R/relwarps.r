#' calculate relative Warp analysis
#' 
#' After Procrustes registration the data is scaled by the bending energy or
#' its inverse to emphasize global/local differences when exploring a sample's
#' shape.
#' 
#' 
#' @param data Input k x m x n real array, where k is the number of points, m
#' is the number of dimensions, and n is the sample size.
#' @param scale Logical: indicating if scaling is requested
#' @param CSinit Logical: if TRUE, all configurations are initially scaled to
#' Unit Centroid Size.
#' @param alpha integer: power of the bending energy matrix. If alpha = 0 then
#' standard Procrustes PCA is carried out. If alpha = 1 then large scale
#' differences are emphasized, if alpha = -1 then small scale variations are
#' emphasised.
#' @param tol tolerance for the eigenvalues of the bending energy matrix to be
#' zero
#' @param orp logical: request orthogonal projection into tangent space.
#' @param pcAlign logical: if TRUE, the shapes are aligned by the principal axis of the first specimen
#' @return
#' \item{bescores }{relative warp scores (PC-scores if \code{alpha = 0})}
#' \item{uniscores }{uniform scores, NULL if  \code{alpha = 0}}
#' \item{Var }{non-affine variation explained by each relative warp}
#' \item{mshape }{sample's conensus shape}
#' \item{rotated }{Procrustes superimposed data}
#' \item{bePCs }{vector basis of nonaffine shape variation- relative warps (plain PCs if  \code{alpha = 0})}
#' \item{uniPCs }{vector basis of affine shape variation - uniform
#' component. NULL if \code{alpha = 0}}
#' @author Stefan Schlager
#' @references Bookstein FL 1989. Principal Warps: Thin-plate splines and the
#' decomposition of deformations. IEEE Transactions on pattern analysis and
#' machine intelligence 11.
#'
#' Bookstein FL, 1991. Morphometric tools for landmark
#' data. Geometry and biology. Cambridge Univ. Press, Cambridge.
#' 
#' Rohlf FJ, Bookstein FL 2003. Computing the Uniform Component of Shape
#' Variation. Systematic Biology 52:66-69.
#' 
#' @examples
#' 
#' data(boneData)
#' pop <- name2factor(boneLM,which=3)
#' rW <- relWarps(boneLM, alpha = -1)
#' \dontrun{
#' require(car)
#' # plot first 5 relative warps scores grouped by population
#' spm(rW$bescores[,1:5],group=pop)
#' # plot uniform component scores grouped by population
#' spm(rW$uniscores[,1:5],group=pop)
#' ##plot non-affine variance associated with each relative warp
#' barplot(rW$Var[,2], xlab="relative Warps")
#' ## visualize first relative warp +-3 sd of the scores
#' rw1 <- showPC(as.matrix(c(-3,3)*sd(rW$bescores[,1])),rW$bePCs[,1,drop=FALSE],rW$mshape)
#' deformGrid3d(rw1[,,1],rw1[,,2],ngrid=5)
#' 
#' ## 2D example:
#' require(shapes)
#' data <- bindArr(gorf.dat, gorm.dat, along=3)
#' sex <- factor(c(rep("fem", dim(gorf.dat)[3]), rep("male",dim(gorm.dat)[3])))
#' rW <- relWarps(data, alpha = -1)
#' # plot first 3 relative warps scores grouped by population
#' spm(rW$bescores[,1:3],group=sex)
#' # plot uniform component scores grouped by population
#' spm(rW$uniscores[,1:2],group=sex)
#' ##plot non-affine variance associated with each relative warp
#' barplot(rW$Var[,2], xlab="relative Warps")
#' ## visualize first relative warp +-3 sd of the scores
#' rw1 <- showPC(as.matrix(c(-3,3)*sd(rW$bescores[,1])),rW$bePCs[,1,drop=FALSE],rW$mshape)
#' deformGrid2d(rw1[,,1],rw1[,,2],ngrid=10)
#' }
#' 
#' @export
relWarps <- function(data,scale=TRUE,CSinit=TRUE,alpha=1,tol=1e-10,orp=TRUE, pcAlign=TRUE)
{
                                        #n <- dim(data)[3]
    uniscores <- uniPCs <- NULL
    m <- dim(data)[2]
    k <- dim(data)[1]
    datanames <- dimnames(data)[[3]]
### superimpose data ###
    proc <- ProcGPA(data,scale=scale,CSinit=CSinit,silent=TRUE,pcAlign=pcAlign)
    if (orp) {
        if (CSinit)
            proc$rotated <- orp(proc$rotated, mshape=proc$mshape)
        else
            message("\n   NOTE: projection into tangent space has been skipped because CSinit == FALSE\n")
    }

    if (alpha !=0 ) {
### create bending energy matrix ###
        BE <- CreateL(proc$mshape,output="Lsubk")$Lsubk
### vectorize and scale superimposed data ###
        vecs <- vecx(proc$rotated)
        vecs <- scale(vecs, scale=FALSE)
        dimnames(proc$rotated)[[3]] <- rownames(vecs) <- datanames
        
### generate covariance matrix of superimposed data ###
        Sc <- cov(vecs)
        
### explore eigenstructure of BE ###	
        eigBE <- svd(BE)
        eigBE$d <- Re(eigBE$d)
        eigBE$v <- Re(eigBE$v)
        
        zero <- which(eigBE$d<tol)
        diaginv <- diagBE <- eigBE$d*0
        diagBE[-zero] <- eigBE$d[-zero]^(-alpha/2)
        diaginv[-zero] <- eigBE$d[-zero]^(alpha/2)
        IM <- Matrix::Diagonal(x=rep(1,m))
        suppressMessages(BE2 <- IM%x%(eigBE$v%*%Matrix::Diagonal(x=diagBE)%*%t(eigBE$v)))
        
### generate covariance structure of scaled space ###
        covcom <- suppressMessages(BE2%*%Sc%*%BE2)
        eigCOVCOM <- svd(covcom)
        nonz <- which(eigCOVCOM$d > tol)
        bescores <- as.matrix(t(suppressMessages(t(eigCOVCOM$v[,nonz])%*%BE2)%*%t(vecs)))[,nonz]
        rownames(bescores) <- rownames(vecs)
        bePCs <-  suppressMessages(IM %x% eigBE$v)
        bePCs <- as.matrix(suppressMessages(bePCs %*% Matrix::Diagonal(x=rep(diaginv,m)) %*% t(bePCs) %*%  eigCOVCOM$v[,nonz]))
       
### calculate uniform component scores ###
        ## U <- NULL
        
        
### Rohlf first method ###
        
        E <- eigBE$v[,-zero]
        N <- diag(rep(1,k))-E%*%solve(crossprod(E))%*%t(E)
        V <- vecs
        NIk <- IM %x% N
        
        svdBend <- svd(V%*%NIk)
        useBendv <- min(ncol(svdBend$v),(m+0.5*m*(m-1)-1))
        LS <- svdBend$u%*%diag(svdBend$d)
        useLS <- min(ncol(LS),(m+0.5*m*(m-1)-1))
        uniscores <- LS[,1:useLS]
        rownames(uniscores) <- datanames
        uniPCs <- svdBend$v[,1:useBendv]
        Var <- createVarTable(eigCOVCOM$d[nonz],square = FALSE)
    } else {
        pca <- prcomp(vecx(proc$rotated))
        bad <- which(pca$sdev^2 < tol)
        bePCs <- pca$rotation[,-bad]
        bescores <- pca$x[,-bad]
        Var <- createVarTable(pca$sdev)
    }
        
### create Variance table according to eigenvalues ###
    
    
    
    
    return(list(bescores=bescores,uniscores=uniscores,Var=Var,mshape=proc$mshape,rotated=proc$rotated,bePCs=bePCs,uniPCs=uniPCs))
    
}

