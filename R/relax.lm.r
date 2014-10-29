#' relax one specific landmark configuration against a reference
#' 
#' relax one specific landmark configuration against a reference (e.g. a
#' sample mean)
#' 
#' 
#' @param lm k x 3 or k x 2 matrix containing landmark data to be slidden.
#' @param reference k x 3 or k x 2 matrix containing landmark of the reference
#' @param SMvector A vector containing the landmarks on the curve(s) that are
#' allowed to slide
#' @param outlines A vector (or if threre are several curves) a list of vectors
#' (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in
#' the successive position on the curve - including the beginning and end
#' points, that are not allowed to slide.
#' @param surp A vector containing Semilandmarks positioned on surfaces.
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
#' @param missing vector of integers, specifying missing (semi-)landmarks. They will be relaxed freely in 3D and not projected onto the target (works only for 2D data).
#' @return returns kx3 matrix of slidden landmarks
#' @author Stefan Schlager
#' @seealso \code{\link{slider3d}}
#' @references Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005.
#' Semilandmarks in Three Dimensions, in Modern Morphometrics in Physical
#' Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer
#' Academic/Plenum Publishers.
#' 
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
#' ## to reduce this example's computation time,
#' # we only use the right hand semi-landmarks
#' # (which keeps the left hand ones fix)
#' surp <- surp[1:316]
#' 
#' relax <- relaxLM(shortnose.lm[1:323, ],
#'          longnose.lm[1:323, ], mesh=shortnose.mesh, iterations=1,
#'          SMvector=fix, deselect=TRUE, surp=surp)
#' 
#' \dontrun{
#' # visualize differences red=before and green=after sliding
#' deformGrid3d(shortnose.lm[1:323, ], relax, ngrid=0)
#' # add surface
#' wire3d(shortnose.mesh, col="white")
#' }
#' 
#' @export
relaxLM <- function(lm,reference,SMvector,outlines=NULL,surp=NULL,sur.name=NULL,mesh=NULL,tol=1e-05,deselect=FALSE,inc.check=TRUE,iterations=0, fixRepro=TRUE, missing=NULL) {
    
    k <- dim(lm)[1]
    m <- dim(lm)[2]
    free <- NULL
    p1 <- 10^12
    lm.orig <- lm
    L <- CreateL(reference)
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
        if (m == 3)
            U <- .calcTang_U_s(vs,vn,SMvector=SMvector,outlines=outlines,surface=surp,deselect=deselect,free=free)
        else
            U <- .calcTang_U(vs,SMvector=SMvector,outlines=outlines,deselect=deselect)
        dataslido <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m)$Gamatrix
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
    return(vs)
}



