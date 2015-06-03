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
#' \item{bescores }{relative warp scores}
#' \item{uniscores }{uniform scores}
#' \item{Var }{non-affine variation explained by each relative warp}
#' \item{mshape }{sample's conensus shape}
#' \item{rotated }{Procrustes superimposed data}
#' \item{bePCs }{vector basis of nonaffine shape variation- relative warps}
#' \item{uniPCs }{vector basis of affine shape variation - uniform
#' component}
#' @author Stefan Schlager
#' @references Bookstein FL 1989. Principal Warps: Thin-plate splines and the
#' decomposition of deformations. IEEE Transactions on pattern analysis and
#' machine intelligence 11. Bookstein FL, 1991. Morphometric tools for landmark
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
#' }
#' 
#' @export
relWarps <- function(data,scale=TRUE,CSinit=TRUE,alpha=1,tol=1e-10,orp=TRUE, pcAlign=TRUE)
{
    #n <- dim(data)[3]
    m <- dim(data)[2]
    k <- dim(data)[1]
    datanames <- dimnames(data)[[3]]
### superimpose data ###
    proc <- ProcGPA(data,scale=scale,CSinit=CSinit,silent=TRUE,pcAlign=pcAlign)
    if (orp)
        proc$rotated <- orp(proc$rotated, mshape=proc$mshape)
    
### create bending energy matrix ###
    BE <- CreateL(proc$mshape,output="Lsubk")$Lsubk
### vectorize and scale superimposed data ###
    vecs <- vecx(proc$rotated)
    vecs <- scale(vecs, scale=FALSE)
    dimnames(proc$rotated)[[3]] <- rownames(vecs) <- datanames
    
### generate covariance matrix of superimposed data ###
    Sc <- cov(vecs)
    
### explore eigenstructure of BE ###	
    eigBE <- eigen(BE,symmetric=TRUE)
    eigBE$values <- Re(eigBE$values)
    eigBE$vectors <- Re(eigBE$vectors)
    
    zero <- which(eigBE$values<tol)
    diaginv <- diagBE <- eigBE$values*0
    diagBE[-zero] <- eigBE$values[-zero]^(-alpha/2)
    diaginv[-zero] <- eigBE$values[-zero]^(alpha/2)
    IM <- diag(rep(1,m))
    
    if (alpha !=0) {	
        BE2 <- IM%x%(eigBE$vectors%*%diag(diagBE)%*%t(eigBE$vectors))
        ##  invBE2 <- IM%x%(eigBE$vectors%*%diag(diaginv)%*%t(eigBE$vectors))
    } else {
        BE2 <- diag(rep(1,k*m))
        ##   invBE2 <- BE2
    }
    
### generate covariance structure of scaled space ###
    covcom <- BE2%*%Sc%*%BE2	
    eigCOVCOM <- eigen(covcom,symmetric=TRUE)
    nonz <- which(eigCOVCOM$values>tol)
    bescores <- t(t(eigCOVCOM$vectors[,nonz])%*%BE2%*%t(vecs))
    rownames(bescores) <- rownames(vecs)
    bePCs <- IM %x% eigBE$vectors
    bePCs <- bePCs %*% diag(rep(diaginv,m)) %*% t(bePCs) %*%  eigCOVCOM$vectors
    
### calculate uniform component scores ###
    ## U <- NULL
    uniscores <- NULL
    
### Rohlf first method ###
    
    E <- eigBE$vectors[,-zero]
    N <- diag(rep(1,k))-E%*%solve(crossprod(E))%*%t(E)
    V <- vecs
    NIk <- IM %x% N
    
    svdBend <- svd(V%*%NIk)
    LS <- svdBend$u%*%diag(svdBend$d)
    uniscores <- LS[,1:(m+0.5*m*(m-1)-1)]
    rownames(uniscores) <- datanames   
    
### create Variance table according to eigenvalues ###
    values <- eigCOVCOM$values[nonz]
    if (length(values)==1) {
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
    
    return(list(bescores=bescores,uniscores=uniscores,Var=Var,mshape=proc$mshape,rotated=proc$rotated,bePCs=bePCs,uniPCs=svdBend$v[,1:(m+0.5*m*(m-1)-1)]))
    
}

