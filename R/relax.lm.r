#' relax one specific landmark configuration against a reference
#' 
#' relax one specific landmark configuration against a reference (e.g. a
#' sample mean)
#' 
#' 
#' @param lm k x 3 or k x 2 matrix containing landmark data to be slidden - or a triangular mesh of class "mesh3d". See details
#' @param reference k x 3 or k x 2 matrix containing landmark of the reference, or a mesh with the same amount of vertices as there are landmarks in \code{lm}.
#' @param SMvector A vector containing the row indices of (semi-) landmarks on the curve(s) that are
#' allowed to slide
#' @param outlines A vector (or if threre are several curves) a list of vectors
#' (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in
#' the successive position on the curve - including the beginning and end
#' points, that are not allowed to slide.
#' @param surp integer vector containing the row indices of semi-landmarks positioned on surfaces.
#' @param sur.name character: containing the filename of the corresponding
#' surface.When specified, mesh has to be NULL.
#' @param mesh triangular mesh of class "mesh3d" loaded into the R workspace,
#' when specified, "sur.name" has to be NULL. The function
#' \code{\link{closemeshKD}} will be used for reprojection onto the surface.
#' @param tol numeric: Threshold for convergence in the sliding proces. Full
#' Procrustes distance between actual result and previous iteration.
#' @param deselect Logical: if TRUE, the SMvector is interpreted as those
#' landmarks, that are not allowed to slide.
#' @param inc.check Logical: if TRUE, the program stops when convergence
#' criterion starts increasing and reports result from last iteration.
#' @param iterations integer: maximum amounts the algorithm runs - even when
#' 'tol' is not reached. When iterations=0, the algorithm runs until
#' convergence.
#' @param fixRepro logical: if \code{TRUE}, fix landmarks will also be
#' projected onto the surface. If you have landmarks not on the surface, select
#' \code{fixRepro=FALSE}
#' @param missing vector of integers, specifying row indices of missing (semi-)landmarks. They will be relaxed freely in 3D and not projected onto the target (works only for 2D data).
#' @param bending if TRUE, bending energy will be minimized, Procrustes distance otherwise (not suggested with large shape differences)
#' @param stepsize integer: dampening factor for the amount of sliding.
#' Useful to keep semi-landmarks from sliding too far off the surface.
#' The displacement is calculated as  \eqn{\Upsilon = \Upsilon^0 + stepsize * UT}{Y = Y0 + stepsize * UT}.
#' Default is set to 1 for bending=TRUE and 0.5 for bending=FALSE.
#' @param use.lm indices specifying a subset of (semi-)landmarks to be used in the rotation step - only used if \code{bending=FALSE}.
#' @param ... additonal arguments - currently unused
#' @return returns kx3 matrix of slidden landmarks
#' @author Stefan Schlager
#' @seealso \code{\link{slider3d}}
#' @references Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005.
#' Semilandmarks in Three Dimensions, in Modern Morphometrics in Physical
#' Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer
#' Academic/Plenum Publishers.
#' @details if \code{lm} is a surface mesh, all vertices will be treated as semilandmarks and a allowed to freely slide along the surface.
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' ### relax shornose against longnose
#' 
#' # define fix landmarks
#' fix <- c(1:5,20:21)
#' # define surface patch by specifying row indices of matrices
#' # all except those defined as fix
#' surp <- c(1:dim(shortnose.lm)[1])[-fix]
#'  
#' relax <- relaxLM(shortnose.lm,
#'          longnose.lm, mesh=shortnose.mesh, iterations=1,
#'          SMvector=fix, deselect=TRUE, surp=surp)
#'
#' ## example minimizing Procrustes distance when displacement is not
#' ## dampened by stepsize
#' relaxProcD <- relaxLM(shortnose.lm,
#'          longnose.lm, mesh=shortnose.mesh, iterations=1,
#'          SMvector=fix, deselect=TRUE, surp=c(1:623)[-fix],bending=FALSE,stepsize=1)
#' 
#' \dontrun{
#' # visualize differences red=before and green=after sliding
#' deformGrid3d(shortnose.lm, relax, ngrid=0)
#'
#'  
#' # visualize differences minimizing Procrusted distances red=before and green=after sliding
#'
#' deformGrid3d(shortnose.lm, relaxProcD, ngrid=0)
#' ## no smooth displacement, now let's check the distances:
#' rot2ref <- rotonto(relaxProcD,longnose.lm)
#' angle.calc(rot2ref$X,rot2ref$Y)
#' # 0.2492027 Procrustes distance between reference and slided shape
#' # (minimizing Procrustes distance)
#' rot2refBend <- rotonto(relax,longnose.lm)
#' angle.calc(rot2refBend$X,rot2refBend$Y)
#' # 0.2861322 Procrustes distance between reference and slided shape
#' # (minimizing bending energy)
#' 
#' rot2refOrig <- rotonto(shortnose.lm,longnose.lm)
#' angle.calc(rot2refOrig$X,rot2refOrig$Y)
#' # 0.3014957 Procrustes distance between reference and original shape
#' ##result: while minimizing Procrustes distance, displacement is not
#' ##guaranteed to be smooth
#' 
#' # add surface
#' wire3d(shortnose.mesh, col="white")
#'
#' 
#' ## finally relax two meshes with corresponding vertices:
#' 
#' mediumnose.mesh <- tps3d(shortnose.mesh,shortnose.lm, (shortnose.lm+longnose.lm)/2)
#' ## we use Procrustes distance as criterion as bending energy is pretty slow because
#' ## of too many coordinates (more than 3000 is very unreasonable).
#' relaxMesh <- relaxLM(shortnose.mesh,mediumnose.mesh,iterations=2,bending=FALSE,stepsize=0.05)
#' }
#' @rdname relaxLM
#' @export
relaxLM <- function(lm,...)UseMethod("relaxLM")

#' @rdname relaxLM
#' @export
relaxLM.matrix <- function(lm,reference,SMvector,outlines=NULL,surp=NULL,sur.name=NULL,mesh=NULL,tol=1e-05,deselect=FALSE,inc.check=TRUE,iterations=0, fixRepro=TRUE, missing=NULL, bending=TRUE,stepsize=ifelse(bending,1,0.5),use.lm=NULL,...) {
    if(inherits(reference,"mesh3d"))
       reference <- vert2points(reference)
    k <- dim(lm)[1]
    m <- dim(lm)[2]
    weights <- NULL
    if (!is.null(use.lm)) {
        weights <- rep(0,nrow(lm))
        weights[use.lm] <- 1
    }
    
    free <- NULL
    p1 <- 10^12
    lm.orig <- lm
    reference <- scale(reference, scale=FALSE)
    if (bending)
        L <- CreateL(reference,output="Lsubk3")

    if (deselect)
        fixLM <- SMvector
    else if (length(SMvector) < k)
        fixLM <- c(1:k)[-SMvector]
    else
        fixRepro <- TRUE

    if (iterations == 0)
        iterations <- 1e10
    if (m == 3) {
        cat(paste("Points will be initially projected onto surfaces","\n","-------------------------------------------","\n"))
        
        if (is.null(mesh)) {
            tmp <- projRead(lm, sur.name)
            
        } else {
            tmp <- projRead(lm,mesh)
        }
        vs <- vert2points(tmp)
        vn <- t(tmp$normals[1:3,])
        
        if (!fixRepro)# use original positions for fix landmarks
            vs[fixLM,] <- lm.orig[fixLM,]
        if (length(missing)) {
            free <- missing
            vs[missing,] <- lm.orig[missing,]
        }
    } else {
        vs <- lm
    }
    count <- 1
    while (p1 > tol && count <= iterations) {
        lm_old <- vs
        cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated
        if (!bending) {
            print(length(weights))
            rot <- rotonto(reference,vs,reflection=FALSE,scale=TRUE,weights=weights,centerweight = TRUE)
            vs <- rot$yrot
            if (m == 3)
                vn <- vn%*%rot$gamm
        }
        if (m == 3)
            U <- .calcTang_U_s(vs,vn,SMvector=SMvector,outlines=outlines,surface=surp,deselect=deselect,free=free)
        else
            U <- .calcTang_U(vs,SMvector=SMvector,outlines=outlines,deselect=deselect)
        if (bending)
            dataslido <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m,stepsize=stepsize)
        else {
            dataslido <- calcProcDGamma(U$U,U$Gamma0,reference,dims=m,stepsize=stepsize)
            dataslido <- rotreverse(dataslido,rot)
        }
        if (m == 3) {
            if (is.null(mesh)) {
                tmp <- projRead(dataslido, sur.name)
            } else {
                tmp <- projRead(dataslido,mesh)
            }
            vs <- vert2points(tmp)
            vn <- t(tmp$normals[1:3,])
            
            if (!fixRepro)# use original positions for fix landmarks
                vs[fixLM,] <- lm.orig[fixLM,]
            if (length(missing))
                vs[missing,] <- dataslido[missing,]
        } else {
            vs <- dataslido
        }
        
        p1_old <- p1
        testproc <- rotonto(lm_old,vs)			   	
        p1 <- sum(diag(crossprod((testproc$X/cSize(testproc$X))-(testproc$Y/cSize(testproc$Y)))))### check for increasing convergence criterion ###		
        if (inc.check) {
            if (p1 > p1_old) {
                vs <- lm_old
                cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
                p1 <- 0
                count <- count+1   
            } else {
                cat(paste("squared distance between iterations:",p1,sep=" "),"\n","-------------------------------------------","\n")
                count <- count+1
            }
        } else {
            cat(paste("squared distance between iterations:",p1,sep=" "),"\n","-------------------------------------------","\n")
            count <- count+1
        }
    }
    gc()
    return(vs)
}

#' @rdname relaxLM
#' @export
relaxLM.mesh3d <- function(lm,reference,tol=1e-05,deselect=FALSE,inc.check=TRUE,iterations=0, fixRepro=TRUE, missing=NULL, bending=FALSE,stepsize=ifelse(bending,1,0.5),use.lm=NULL, ...){
    lmtmp <- vert2points(lm)
    mesh <- lm
    SMvector <- surp <- 1:ncol(lm$vb)
    sur.name <- NULL
    outlines <- NULL
    reltmp <- relaxLM(lmtmp,reference=reference,SMvector=SMvector,surp=surp,mesh=mesh,tol=tol,deselect=deselect,inc.check=inc.check,iterations=iterations, fixRepro=fixRepro, missing=missing, bending=bending,stepsize=stepsize,use.lm=use.lm)
    lm$vb[1:3,] <- t(reltmp)
    lm <- vcgUpdateNormals(lm)
    return(lm)
}
