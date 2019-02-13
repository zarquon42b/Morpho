#' slides Semilandmarks along curves and surfaces in 3D by minimising bending
#' energy of a thin-plate spline deformation.
#' 
#' slides Semilandmarks along curves and surfaces in 3D. The positions on the
#' surface are sought which minimise bending energy (of a thin-plate spline
#' deformation)
#' 
#' 
#' @param dat.array Input k x m x n real array, where k is the number of
#' points, m is the number of dimensions, and n is the sample size. Ideally the
#' dimnames[[3]] vector contains the names of the surface model (without file
#' extension) - e.g. if the model is named "surface.ply", the name of the
#' corresponding matrix of the array would be "surface"
#' @param SMvector A vector containing the row indices of (semi-) landmarks on the curve(s) and
#' surfaces that are allowed to slide
#' @param outlines A vector (or if threre are several curves) a list of vectors
#' (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in
#' the successive position on the curve - including the beginning and end
#' points, that are not allowed to slide.
#' @param surp integer vector containing the row indices of semi-landmarks positioned on surfaces.
#' @param sur.path Path to the surface models (e.g. ply, obj, stl files)
#' @param sur.name character vector: containing the filenames of the
#' corresponding surfaces - e.g. if the dat.array[,,i] belongs to
#' surface_i.ply, sur.name[i] would be surface_i.ply. Only necessary if
#' dat.array does not contain surface names.
#' @param meshlist list containing triangular meshes of class 'mesh3d', for
#' example imported with \code{\link{mesh2ply}} or \code{\link{file2mesh}} in
#' the same order as the specimen in the array (see examples below).
#' @param ignore vector containing indices of landmarks that are to be ignored.
#' Indices of outlines/surfaces etc will be updated automatically.
#' @param sur.type character:if all surfaces are of the same file format and
#' the names stored in dat.array, the file format will be specified here.
#' @param tol numeric: Threshold for convergence in the sliding process
#' @param deselect Logical: if TRUE, the SMvector is interpreted as those
#' landmarks, that are not allowed to slide.
#' @param inc.check Logical: if TRUE, the program stops when convergence
#' criterion starts increasing and reports result from last iteration.
#' @param fullGPA Logical: if FALSE, only a partial procrustes fit will be
#' performed.
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
#' @param mc.cores integer: determines how many cores to use for the
#' computation. The default is autodetect. But in case, it doesn't work as
#' expected cores can be set manually.
#' @param fixRepro logical: if \code{TRUE}, fix landmarks will also be
#' projected onto the surface. If you have landmarks not on the surface, select
#' \code{fixRepro=FALSE}
#' @param missingList a list of length samplesize containing integer vectors of row indices specifying missing landmars for each specimen. For specimens without missing landmarks enter \code{numeric(0)}.
#' @param use.lm indices specifying a subset of (semi-)landmarks to be used in the rotation step - only used if \code{bending=FALSE}.
#' @param silent logical: if TRUE, console output is suppressed.
#' @return
#' \item{dataslide }{array containing slidden Landmarks in the original
#' space - not yet processed by a Procrustes analysis}
#' \item{vn.array }{array containing landmark normals}
#' @note if \code{sur.path = NULL} and \code{meshlist = NULL}, surface landmarks
#' are relaxed based on a surface normals approximated by the pointcloud, this can lead to bad results for sparse sets of semilandmarks. Obviously, no projection onto the surfaces will be occur and landmarks will likely be off
#' the original surface.
#' @section Warning: Depending on the size of the suface meshes and especially
#' the amount of landmarks this can use an extensive amount of your PC's
#' resources, especially when running in parallel. As the computation time and
#' RAM usage of matrix algebra involved is quadratic to the amount of landmarks
#' used, doubling the amount of semi-landmarks will quadruple computation time
#' and system resource usage. You can easily stall you computer with this
#' function with inappropriate data.
#' @author Stefan Schlager
#' @seealso \code{\link{relaxLM}, \link{createMissingList}}
#' @encoding utf8
#' @references Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape analysis
#' of symmetric structures: quantifying variation among individuals and
#' asymmetry. Evolution 56(10):1909-1920.
#' 
#' Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. Semilandmarks in Three
#' Dimensions, in Modern Morphometrics in Physical Anthropology. Edited by D.
#' E. Slice, pp. 73-98. New York: Kluwer Academic/Plenum Publishers.
#' 
#' Schlager S. 2012. Sliding semi-landmarks on symmetric structures in three
#' dimensions. American Journal of Physical Anthropology, 147(S52):261. URL:
#' http://dx.doi.org/10.1002/ajpa.21502.
#' 
#' Schlager S. 2013. Soft-tissue reconstruction of the human nose: population
#' differences and sexual dimorphism. PhD thesis,
#' \enc{Universitätsbibliothek}{Universitaetsbibliothek} Freiburg.  URL:
#' \url{http://www.freidok.uni-freiburg.de/volltexte/9181/}.
#' @examples
#' \dontrun{
#' data(nose)
#' ###create mesh for longnose
#' longnose.mesh <- tps3d(shortnose.mesh,shortnose.lm,longnose.lm,threads=1)
#' ### write meshes to disk
#' mesh2ply(shortnose.mesh, filename="shortnose")
#' mesh2ply(longnose.mesh, filename="longnose")
#' 
#' ## create landmark array
#' data <- bindArr(shortnose.lm, longnose.lm, along=3)
#' dimnames(data)[[3]] <- c("shortnose", "longnose")
#' 
#' # define fix landmarks
#' fix <- c(1:5,20:21)
#' # define surface patch by specifying row indices of matrices
#' # all except those defined as fix
#' surp <- c(1:nrow(shortnose.lm))[-fix]
#' 
#' slide <- slider3d(data, SMvector=fix, deselect=TRUE, surp=surp,
#'                   sur.path=".",iterations=1,mc.cores=1)
#'                   # sur.path="." is the current working directory
#' 
#' # now one example with meshes in workspace
#' 
#' meshlist <- list(shortnose.mesh,longnose.mesh)
#' 
#' slide <- slider3d(data, SMvector=fix, deselect=TRUE, surp=surp,
#'                   iterations=1, meshlist=meshlist,
#'                   mc.cores=1,fixRepro=FALSE)
#' require(rgl)
#' ## visualize sliding
#' deformGrid3d(slide$dataslide[,,1],shortnose.lm,ngrid = 0)
#' ## these are fix
#' spheres3d(slide$dataslide[fix,,1],col=4,radius=0.7)
#'
#' ###finally an example with missing landmarks:
#' ## we assume that coordinates 185:189, 205:209 and 225:229 are in the second config are missing
#' missingList <- createMissingList(2)
#' missingList[[2]] <- c(185:189,205:209,225:229)
#' slideMissing <- slider3d(data, SMvector=fix, deselect=TRUE, surp=surp,
#'                   iterations=1, meshlist=meshlist,
#'                   mc.cores=1,fixRepro=FALSE,missingList=missingList)
#'
#' ## example with two curves
#' ## Example with surface semilandmarks and two curves
#' fix <- c(1:5,20:21)
#' outline1 <- c(304:323)
#' outline2 <- c(604:623)
#' outlines <- list(outline1,outline2)
#' surp <- c(1:623)[-c(fix,outline1,outline2)]
#' slideWithCurves <- slider3d(data, SMvector=fix, deselect=TRUE, surp=surp,
#'                             meshlist=meshlist,iterations=1,mc.cores=1,outlines=outlines)
#' deformGrid3d(slideWithCurves$dataslide[,,1],shortnose.lm,ngrid = 0)
#' plot(slideWithCurves)
#'
#' ## finally an example with sliding without meshes by estimating the surface from the
#' ## semi-landmarks
#'
#' slideWithCurvesNoMeshes <- slider3d(data, SMvector=fix, deselect=TRUE, surp=surp,
#'                             iterations=1,mc.cores=1,outlines=outlines)
#' ## compare it to the data with surfaces
#' deformGrid3d(slideWithCurves$dataslide[,,1],slideWithCurvesNoMeshes$dataslide[,,1],ngrid = 0)
#' ## not too bad, only lonely surface semi-landmarks are a bit off
#'
#' }
#' 
#' @export
slider3d <- function(dat.array,SMvector,outlines=NULL,surp=NULL,sur.path=NULL,sur.name=NULL, meshlist=NULL, ignore=NULL,sur.type="ply",tol=1e-05,deselect=FALSE,inc.check=TRUE,recursive=TRUE,iterations=0,initproc=TRUE,fullGPA=FALSE,pairedLM=0,bending=TRUE,stepsize=ifelse(bending,1,0.5),mc.cores = parallel::detectCores(), fixRepro=TRUE,missingList=NULL,use.lm=NULL,silent=FALSE)
{
    if (.Platform$OS.type == "windows" && mc.cores > 1) {
        cl <- makeCluster(mc.cores)            
        registerDoParallel(cl=cl)
    } else if (mc.cores > 1) {
        registerDoParallel(cores = mc.cores)
    } else
        registerDoSEQ()
    
    if (iterations == 0)
        iterations <- 1e10
    
    if (is.null(outlines) && is.null(surp))	
        stop("nothing to slide")
    fixLM <- integer(0)
    n <- dim(dat.array)[3]
    k <- dim(dat.array)[1]
    m <- dim(dat.array)[2]
    nomesh <- FALSE
    if (length(sur.name) && is.null(sur.path))
        sur.path <- ""
    if (is.null(meshlist) && is.null(sur.path))
        nomesh <- TRUE
    if (pairedLM[1]!=0 && is.vector(pairedLM))# check if there are only 2 symmetric lms
        pairedLM <- t(as.matrix(pairedLM))
    if(!is.null(missingList))
        if(length(missingList) != n)
            stop(paste0("missingList must be of length", n," - same as samplesize"))
### update indexing for after ignored landmarks are removed ###	
    if (!is.null(ignore)) {
        li <- length(ignore)
        lm.old <- c(1:k)[-ignore]
        mat.ptr <- matrix(c(1:(k-li),lm.old),k-li,2)
        ptr <- function(xo)	### define pointer function for indexing
        {
            if (length(which(ignore %in% xo))!= 0)
                xo <- xo[-which(xo %in% ignore)]
            for (i in 1:(k-li))
                xo[which(xo==mat.ptr[i,2])] <- mat.ptr[i,1]
            return(xo)
        }
        if (!is.null(missingList))
            missingList <- lapply(missingList,ptr)
        if (!is.null(outlines)) ### update outline indices
            outlines <- lapply(outlines,ptr)
        if (!is.null(surp)) 	### update surface indices
            surp <- ptr(surp)
        
        if (!is.null(SMvector)) ### of fixed/sliding definition
            SMvector <- ptr(SMvector)
        
        if (pairedLM[1]!=0){	### update paired landmarks indices
            count <- 0
            del <- NULL
            for (i in 1:dim(pairedLM)[1]) {	
                if (length(which(ignore %in% pairedLM[i,]))!=0) {
                    count <- count+1
                    del[count] <- i
                }
            }
            pairedLM <- pairedLM[-del,]
            if (is.vector(pairedLM))
                pairedLM <- t(as.matrix(pairedLM))
            
            if (dim(pairedLM)[1]==0) {
                pairedLM <- 0
            } else {
                pairedLM <- apply(pairedLM,2,ptr)
                if (is.vector(pairedLM))
                    pairedLM <- t(as.matrix(pairedLM))
            }
        }
        dat.array <- dat.array[-ignore,,]
        k <- dim(dat.array)[1]
    }
    
    vn.array <- dat.array
    data.orig <- dat.array
    if (deselect)
        fixLM <- SMvector
    else if (length(SMvector) < k)
        fixLM <- c(1:k)[-SMvector]
    else
        fixRepro <- TRUE

    weights <- NULL
    if (!is.null(use.lm)) {
        weights <- rep(0,dim(dat.array)[1])
        weights[use.lm] <- 1
    }
    if(length(sur.name)==0 && !is.null(sur.path)) {
        sur.name <- dimnames(dat.array)[[3]]
        sur.name <- paste(sur.path,"/",sur.name,".",sur.type,sep="")
        
    }
    p1 <- 10^12

    ## check for existing surfaces (if needed)
    if (length(sur.name)) {
        if (length(sur.name) != dim(dat.array)[3])
            stop("length of sur.name does not match number of landmark configurations")
        checkfiles <- file.exists(sur.name)
        if (prod(checkfiles) == 0) {
            warning(paste0("missing mesh files: ", sur.name[checkfiles == 0]))
            stop("Some mesh files do not exist, please check your data and/or correct array naming")
        }
    }
    if (length(meshlist)) {
        if (length(meshlist) != dim(dat.array)[3])
            stop("length of meshlist does not match number of landmark configurations")
    }
    
    ini <- rotonto(dat.array[,,1],dat.array[,,2],signref=FALSE) # create mean between first tow configs to avoid singular BE Matrix
    mshape <- (ini$Y+ini$X)/2
    if (!silent && !nomesh)
        cat(paste("Points will be initially projected onto surfaces","\n","-------------------------------------------","\n"))
    ## parallel function in case meshlist != NULL
    parfunmeshlist <- function(i,data) {
        if (!is.list(data))
            tmpdata <- data[,,i]
        else
            tmpdata <- data[[i]]
        out <- projRead(tmpdata,meshlist[[i]])
        if (!is.null(missingList))
            if(length(missingList[[i]]))
                out$vb[1:3,missingList[[i]]] <- t(tmpdata[missingList[[i]],])
        return(out)
    }
    parfunmeshfile <- function(i, data) {
        if (!is.list(data))
            tmpdata <- data[,,i]
        else
            tmpdata <- data[[i]]
        
        out <- projRead(tmpdata,sur.name[i])
        if (!is.null(missingList))
            if(length(missingList[[i]]))
                out$vb[1:3,missingList[[i]]] <- t(tmpdata[missingList[[i]],])
        
        return(out)

    }
    parfunnomesh <- function(i, data) {
        if (!is.list(data))
            tmpdata <- data[,,i]
        else
            tmpdata <- data[[i]]
        
        out <- vcgUpdateNormals(tmpdata,silent=TRUE)
        return(out)
    }
    
    if (is.null(meshlist) && !nomesh) {
        repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshfile(i,dat.array)
    } else if (!nomesh) {
        repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshlist(i,dat.array)
        
    } else {
        repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunnomesh(i,dat.array)
        message("no surfaces specified - surface is approximated from point cloud")
    }    
    for (j in 1:n) {
        reprotmp <- repro[[j]]         
        dat.array[,,j] <- t(reprotmp$vb[1:3,])
        vn.array[,,j] <- t(reprotmp$normals[1:3,])
        
    }
    
    
    if (!fixRepro)# use original positions for fix landmarks
        dat.array[fixLM,,] <- data.orig[fixLM,,]
    
    if (!silent && !nomesh)
        cat(paste("\n","-------------------------------------------","\n"),"Projection finished","\n","-------------------------------------------","\n")
    
    if (initproc==TRUE) { # perform proc fit before sliding
        if (!silent)
            cat("Inital procrustes fit ...")	
        procini <- ProcGPA(dat.array,scale=fullGPA,silent=silent)
        mshape <- procini$mshape
    }
    dataslide <- dat.array
    
    if (pairedLM[1]!=0) {# create symmetric mean to get rid of assymetry along outlines/surfaces after first relaxation
        Mir <- diag(c(-1,1,1))
        A <- mshape
        Amir <- mshape%*%Mir
        Amir[c(pairedLM),] <- Amir[c(pairedLM[,2:1]),]
        symproc <- rotonto(A,Amir)
        mshape <- (symproc$X+symproc$Y)/2
    }
    if (!silent)
        cat(paste("Start sliding...","\n","-------------------------------------------","\n"))
    gc(verbose=F)
    ## calculation for a defined max. number of iterations
    count <- 1
    while (p1>tol && count <= iterations) {
        dataslide_old <- dataslide
        mshape_old <- mshape
        if (!silent)
            cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated
        if (recursive==TRUE)    # slided Semilandmarks are used in next iteration step
            dat.array <- dataslide
        if (bending)
            L <- CreateL(mshape,output="Lsubk3")
        else
            fixRepro=TRUE
        a.list <- as.list(1:n)
        slido <- function(j)          		
        {
            free <- NULL
            if (!is.null(missingList))
                if(length(missingList[[j]])) {
                    if (!is.null(outlines)) {
                        notoutlines <- which(!missingList[[j]] %in% unlist(outlines))
                        if (length(notoutlines))
                            free <- missingList[[j]][notoutlines]
                    } else {
                        free <- missingList[[j]]
                    }
                }
            tmpdata <- dat.array[,,j]
            tmpvn <- vn.array[,,j]
            if (!bending) {
                rot <- rotonto(mshape,tmpdata,reflection=FALSE,scale=TRUE,weights=weights,centerweight=TRUE)
                tmpdata <- rot$yrot
                tmpvn <- tmpvn%*%rot$gamm
            }
            U <- .calcTang_U_s(tmpdata,tmpvn,SMvector=SMvector,outlines=outlines,surface=surp,deselect=deselect,free=free)
            if (bending) {
                dataslido <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m,stepsize=stepsize)
            } else {
                dataslido <- calcProcDGamma(U$U,U$Gamma0,mshape,dims=m,stepsize=stepsize)
                dataslido <- rotreverse(dataslido,rot)
            }
            return(dataslido)
        }
        
        a.list <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.export=c("calcGamma",".calcTang_U_s"),.packages=c("Morpho","Rvcg")) %dopar% slido(i)
        
        
###projection onto surface
        
        if (is.null(meshlist) && !nomesh) {
            repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshfile(i,a.list)  
        } else if (!nomesh) {
            repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshlist(i,a.list)
        } else {
            repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunnomesh(i,a.list)
        }
        
        for (j in 1:n) {
            reprotmp <- repro[[j]]         
            dataslide[,,j] <- t(reprotmp$vb[1:3,])
            vn.array[,,j] <- t(reprotmp$normals[1:3,])
        }
        
        if (!fixRepro)# use original positions for fix landmarks
            dataslide[fixLM,,] <- data.orig[fixLM,,]

        if (!silent)
            cat("estimating sample mean shape...")          	
        proc <- ProcGPA(dataslide,scale=fullGPA,silent=silent)
        mshape <- proc$mshape
        if (pairedLM[1]!=0) {# create symmetric mean to get rid of assymetry along outline after first relaxation
            Mir <- diag(c(-1,1,1))
            A <- mshape
            Amir <- mshape%*%Mir
            Amir[c(pairedLM),] <- Amir[c(pairedLM[,2:1]),]
            symproc <- rotonto(A,Amir)
            mshape <- (symproc$X+symproc$Y)/2
        }     
        p1_old <- p1
        testproc <- rotonto(mshape_old,mshape)			   	
        p1 <- sum(diag(crossprod((testproc$X/cSize(testproc$X))-(testproc$Y/cSize(testproc$Y)))))
        
### check for increasing convergence criterion ###		
        if (inc.check) {
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
        } else {
            if (!silent)
                cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
            count <- count+1         
        }
        gc(verbose = FALSE)
    }
    gc(verbose = FALSE)
    out <- list(dataslide=dataslide,vn.array=vn.array)
    class(out) <- "slider3d"
    out$sliderinfo <- list(fixLM=fixLM,outlineLM=outlines,surfaceLM=surp)
    if (.Platform$OS.type == "windows" && mc.cores > 1)
        stopCluster(cl)
    return(out)
}
#' plot the result of slider3d
#'
#' plot the result of slider3d
#' @param x result of \code{slider3d} call
#' @param cols vector containing colors for each coordinate type cols[1]=landmarks, cols[2]=surface landmarks, cols[3]=outlines.
#' @param pt.size size of plotted points/spheres. If \code{point="s"}.
#' \code{pt.size} defines the radius of the spheres. If \code{point="p"} it
#' sets the variable \code{size} used in \code{point3d}.
#' @param point how to render landmarks.
#' @param specimen integer: select the specimen to plot
#' @param add logical: if TRUE, a new rgl window is opened.
#' @param ... additonal, currently unused parameters
#' @method plot slider3d
#' @export
plot.slider3d <- function(x, cols=2:4, pt.size = NULL, point = c("sphere", "point"),specimen=1,add = TRUE,...) {
    if (!inherits(x, "slider3d")) 
        stop("please provide object of class 'slider3d'")
    sliderinfo <- x$sliderinfo
    plotx <- x$dataslide[,,specimen]
    point <- match.arg(point[1],c("sphere", "point"))
    radius <- pt.size
    if (is.null(radius)) {
        if (point == "sphere") 
            radius <- (cSize(plotx)/sqrt(nrow(plotx))) * (1/30)
        else radius <- 10
    }
    size <- radius
    if (point == "sphere") 
        rendpoint <- spheres3d
    else if (point == "point") 
        rendpoint <- points3d
    
    if (!add) 
        open3d()
    if (!is.null(sliderinfo$fixLM))
        rendpoint(plotx[sliderinfo$fixLM,], col=cols[1], radius = radius, size = size)
    if (!is.null(sliderinfo$surfaceLM))
        rendpoint(plotx[sliderinfo$surfaceLM,], col=cols[2], radius = radius/2, size = size/2)
    if (!is.null(sliderinfo$outlineLM))
        rendpoint(plotx[unlist(sliderinfo$outlineLM),], col=cols[3], radius = radius/2, size = size/2)
}
