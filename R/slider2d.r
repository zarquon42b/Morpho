#' slides Semilandmarks along curves 2D by minimising bending
#' energy of a thin-plate spline deformation.
#' 
#' slides Semilandmarks along curves 2D. The positions are sought by minimising bending energy (of a thin-plate spline
#' deformation) or Procrustes distance
#' 
#' 
#' @param dataframe Input k x 2 x n real array, where k is the number of
#' points and n is the sample size. Ideally the
#' @param SMvector A vector containing the row indices of (semi-) landmarks on the curve(s) and
#' surfaces that are allowed to slide
#' @param outlines A vector (or if threre are several curves) a list of vectors
#' (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in
#' the successive position on the curve - including the beginning and end
#' points, that are not allowed to slide.
#' @param tol numeric: Threshold for convergence in the sliding process
#' @param deselect Logical: if TRUE, the SMvector is interpreted as those
#' landmarks, that are not allowed to slide.
#' @param recursive Logical: if TRUE, during the iterations of the sliding
#' process, the outcome of the previous iteration will be used.  Otherwise the
#' original configuration will be used in all iterations.
#' @param iterations integer: select manually the max. number of iterations
#' that will be performed during the sliding process (usefull, when there is
#' very slow convergence).  0 means iteration until convergence.
#' @param initproc requests initial Procrustes fit before sliding.
#' @param pairedLM A X x 2 numeric matrix with the indices of the rows
#' containing paired Landmarks. E.g. the left column contains the lefthand
#' landmarks, while the right side contains the corresponding right hand
#' landmarks. - This will ideally create symmetric mean to get rid of
#' assymetry.
#' @param bending if TRUE, bending energy will be minimized, Procrustes distance otherwise.
#' @param stepsize integer: dampening factor for the amount of sliding.
#' Useful to keep semi-landmarks from sliding too far off the surface.
#' The displacement is calculated as  \eqn{\Upsilon = \Upsilon^0 + stepsize * UT}{Y = Y0 + stepsize * UT}.
#' Default is set to 1 for bending=TRUE and 0.5 for bending=FALSE.
#' @param silent logical: if TRUE, console output is suppressed.
#' @return
#' returns an array containing slided coorndinates in the original
#' space - not yet processed by a Procrustes analysis.
#' @section Warning: Depending on the amount of landmarks this can use an extensive amount of your PC's
#' resources, especially when running in parallel. As the computation time and
#' RAM usage of matrix algebra involved is quadratic to the amount of landmarks
#' used, doubling the amount of semi-landmarks will quadruple computation time
#' and system resource usage. You can easily stall you computer with this
#' function with inappropriate data.
#' @author Stefan Schlager
#' @seealso \code{\link{relaxLM}, \link{slider3d}}
#' @encoding utf8
#' @export
slider2d <- function(dataframe,SMvector,outlines,tol=1e-05,deselect=FALSE,recursive=TRUE,iterations=0,initproc=FALSE,pairedLM=NULL,bending=TRUE,stepsize=1,silent=FALSE)
{
    n <- dim(dataframe)[3]
    k <- dim(dataframe)[1]
    m <- dim(dataframe)[2]
    p1 <- 10^12
    if (iterations == 0)
        iterations <- 1e10
    
    ini <- rotonto(dataframe[,,1],dataframe[,,2],reflection=T)
    mshape <- (ini$X+ini$Y)/2
    
    if(initproc==TRUE) { # perform proc fit before sliding
        procini <- ProcGPA(dataframe,scale=TRUE,silent=TRUE)
        mshape <- procini$mshape
    }
    dataslide <- dataframe
    
    if (!is.null(pairedLM)) {# create symmetric mean to get rid of assymetry along outline after first relaxation
        Mir <- diag(c(-1,1,1))
        A <- mshape
        Amir <- mshape%*%Mir
        Amir[c(pairedLM),] <- Amir[c(pairedLM[,2:1]),]
        symproc <- rotonto(A,Amir)
        mshape <- (symproc$X+symproc$Y)/2
    }
    count <- 1
    while (p1>tol && count <= iterations) {
        dataslide_old <- dataslide
        mshape_old <- mshape
        if (!silent) 
            cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated  
        
        if (recursive==TRUE)      # slided Semilandmarks are used in next iteration step
            dataframe <- dataslide
        
        L <- CreateL(mshape)
        
        for (j in 1:n) {
            tmp <- dataframe[,,j]
            if (!bending) {
                rot <- rotonto(mshape,tmp,scale=TRUE,reflection=FALSE)
                tmp <- rot$yrot
            }
            U <- .calcTang_U(tmp,SMvector=SMvector,outlines=outlines,deselect=deselect)
            if (bending) {
                dataslide[,,j] <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m,stepsize=stepsize)
            } else {
                tmpslide <- calcProcDGamma(U$U,U$Gamma0,mshape,dims=m,stepsize=stepsize)
                dataslide[,,j] <- rotreverse(tmpslide,rot)
            }
        }
        proc <- ProcGPA(dataslide,scale=TRUE,silent=TRUE)
        mshape <- proc$mshape
        p1_old <- p1   
        p1 <- sum(diag(crossprod((mshape_old/cSize(mshape_old))-(mshape/cSize(mshape)))))
                       
        ## check for increasing convergence criterion ###		
        if (p1 > p1_old) {
            dataslide <- dataslide_old
            if (!silent)
                cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
            p1 <- 0
        } else {
            if (!silent)
                cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
            count <- count+1 
        }          		
    }
    return(dataslide)
}
