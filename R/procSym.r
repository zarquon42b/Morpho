#' Procrustes registration
#' 
#' \code{procSym} performs Procrustes superimposition including sliding of
#' semi-landmarks on curves/outlines in 2D and 3D.
#' 
#' This function performs Procrustes registration, allowing a variety of
#' options, including scaling, orthogonal projection into tangentspace and
#' relaxation of semi-landmarks on curves (without reprojection onto the
#' surface/actual outline). It also allows the superimpositioning to be
#' performed using only a subset of the available landmark.  For taking into
#' account object symmetry, \code{pairedLM} needs to be set. This generates an
#' object of class \code{"symproc"}. Otherwise an object of class
#' \code{"nosymproc"}.
#' 
#' @param dataarray Input k x m x n real array, where k is the number of
#' points, m is the number of dimensions, and n is the sample size.
#' @param scale logical: indicating if scaling is requested to minimize the General Procrustes distance. To avoid all scaling, one has to set \code{CSinit=FALSE}, too.
#' @param reflect logical: allow reflections.
#' @param CSinit logical: if TRUE, all configurations are initially scaled to
#' Unit Centroid Size.
#' @param orp logical: if TRUE, an orthogonal projection at the meanshape into
#' tangent space is performed.
#' @param tol numeric: Threshold for convergence in the sliding process
#' @param pairedLM A X x 2 matrix containing the indices (rownumbers) of the
#' paired LM. E.g. the left column contains the lefthand landmarks, while the
#' right side contains the corresponding right hand landmarks.
#' @param sizeshape Logical: if TRUE, a log transformed variable of Centroid
#' Size will be added to the shapedata as first variable before performing the
#' PCA.
#' @param use.lm vector of integers to define a subset of landmarks to be used
#' for Procrustes registration.
#' @param center.part Logical: if TRUE, the data superimposed by the subset
#' defined by use.lm will be centered according to the centroid of the complete
#' configuration. Otherwise orp will be set to FALSE to avoid erroneous
#' projection into tangent space.
#' @param pcAlign logical: if TRUE, the shapes are aligned by the principal axis of the first specimen
#' @param distfun character: "riemann" requests a Riemannian distance for
#' calculating distances to mean, while "angle" uses an approximation by
#' calculating the angle between rotated shapes on the unit sphere.
#' @param SMvector A vector containing the landmarks on the curve(s) that are
#' allowed to slide
#' @param outlines A vector (or if threre are several curves) a list of vectors
#' (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in
#' the successive position on the curve - including the beginning and end
#' points, that are not allowed to slide.
#' @param deselect Logical: if TRUE, the SMvector is interpreted as those
#' landmarks, that are not allowed to slide.
#' @param recursive Logical: if TRUE, during the iterations of the sliding
#' process, the outcome of the previous iteration will be used.  Otherwise the
#' original configuration will be used in all iterations.
#' @param iterations integer: select manually how many iterations will be
#' performed during the sliding process (usefull, when there is very slow
#' convergence).  0 means iteration until convergence.
#' @param initproc Logical: indicating if the first Relaxation step is
#' performed against the mean of an initial Procrustes superimposition.
#' Symmetric configurations will be relaxed against a perfectly symmetrical mean.
#' @param bending if TRUE, bending energy will be minimized, Procrustes distance otherwise (not suggested with large shape differences)
#' @param stepsize integer: dampening factor for the sliding.
#' Useful to keep semi-landmarks from sliding too far off the surface.
#' The displacement is calculated as \cr
#' \code{stepsize * displacement}.
#' 
#' @return
#' \item{size }{a vector containing the Centroid Size of the configurations}
#' \item{rotated }{k x m x n array of the rotated configurations}
#' \item{Sym }{k x m x n array of the Symmetrical component - only
#' available for the "Symmetry"-Option (when pairedLM is defined)}
#' \item{Asym }{k x m x n array of the Asymmetrical component - only
#' available for the "Symmetry"-Option (when pairedLM is defined)}
#' \item{asymmean }{k x m matrix of mean asymmetric deviation from
#' symmetric mean}
#' \item{mshape }{sample meanshape}
#' \item{symmean }{meanshape of symmetrized configurations}
#' \item{tan }{if orp=TRUE: Residuals in tangentspace else, Procrustes
#' residuals - only available without the "Symmetrie"-Option}
#' \item{PCs }{Principal Components - if sizeshape=TRUE, the first variable
#' of the PCs is size information (as log transformed Centroid Size)}
#' \item{PCsym }{Principal Components of the Symmetrical Component}
#' \item{PCasym }{Principal Components of the Asymmetrical Component}
#' \item{PCscores }{PC scores}
#' \item{PCscore_sym }{PC scores of the Symmetrical Component}
#' \item{PCscore_asym }{PC scores of the Asymmetrical Component}
#' \item{eigenvalues }{eigenvalues of the Covariance matrix}
#' \item{eigensym }{eigenvalues of the "Symmetrical" Covariance matrix}
#' \item{eigenasym }{eigenvalues of the "Asymmetrical" Covariance matrix}
#' \item{Variance }{Table of the explained Variance by the PCs}
#' \item{SymVar }{Table of the explained "Symmetrical" Variance by the PCs}
#' \item{AsymVar }{Table of the explained "Asymmetrical" Variance by the PCs}
#' \item{orpdata }{k x m x n array of the rotated configurations projected
#' into tangent space}
#' \item{rho }{vector of Riemannian distance from the mean}
#' \item{dataslide }{array containing slidden Landmarks in the original
#' space - not yet processed by a Procrustes analysis. Only available if a
#' sliding process was requested}
#' @note For processing of surface landmarks or including the reprojection of
#' slid landmarks back onto 3D-surface representations, the usage of
#' \code{\link{slider3d}} is recommended.
#' @author Stefan Schlager
#' @seealso \code{\link{slider3d}}
#' @references Dryden IL, and Mardia KV. 1998. Statistical shape analysis.
#' Chichester.
#' 
#' Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape analysis of symmetric
#' structures: quantifying variation among individuals and asymmetry. Evolution
#' 56(10):1909-1920.
#' 
#' Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. Semilandmarks in Three
#' Dimensions, in Modern Morphometrics in Physical Anthropology. Edited by D.
#' E. Slice, pp. 73-98. New York: Kluwer Academic/Plenum Publishers.
#' 
#' @examples
#' 
#' require(rgl)
#' data(boneData)
#' 
#' ### do an analysis of symmetric landmarks
#' ## visualize landmarks on surface
#' \dontrun{
#'  spheres3d(boneLM[,,1])
#' wire3d(skull_0144_ch_fe.mesh,col=3)
#' ## get landmark numbers
#' text3d(boneLM[,,1],text=paste(1:10),adj = 1, cex=3)
#' }
#' ## determine paired Landmarks left side:
#' left <- c(4,6,8)
#' ## determine corresponding Landmarks on the right side:
#' # important: keep same order
#' right <- c(3,5,7)
#' pairedLM <- cbind(left,right)
#' symproc <- procSym(boneLM, pairedLM=pairedLM)
#' \dontrun{
#' ## visualize first 3 PCs of symmetric shape
#' pcaplot3d(symproc, sym=TRUE)
#' ## visualize first 3 PCs of asymmetric shape
#' pcaplot3d(symproc, sym=FALSE)
#' 
#' ## visualze distribution of symmetric PCscores population
#' pop <- name2factor(boneLM, which=3)
#' require(car)
#' spm(~symproc$PCscore_sym[,1:5], groups=pop)
#' ## visualze distribution of asymmetric PCscores population
#' spm(~symproc$PCscore_asym[,1:5], groups=pop)
#' }
#' 
#' 
#' @export
procSym <- function(dataarray, scale=TRUE, reflect=TRUE, CSinit=TRUE,  orp=TRUE, tol=1e-05, pairedLM=NULL, sizeshape=FALSE, use.lm=NULL, center.part=FALSE, pcAlign=TRUE, distfun=c("angle", "riemann"), SMvector=NULL, outlines=NULL, deselect=FALSE, recursive=TRUE,iterations=0, initproc=FALSE, bending=TRUE,stepsize=1)
{
    t0 <- Sys.time()     
    A <- dataarray
    k <- dim(A)[1]
    m <- dim(A)[2]     
    n <- dim(A)[3]
    Mir <- rep(1, m)
    Mir[1] <- -1
    Mir <- diag(Mir)
    dataslide <- NULL
    CS <- NULL
    if (substr(distfun[1], 1L, 1L) == "r"){
        distfun <- kendalldist
    } else {
        distfun <- function(x,y)
            {
                rho <- angle.calc(x,y)
                return(rho)
            }
    }
    if (is.null(SMvector)) 
        CS <- apply(A, 3, cSize)
    if (!is.null(SMvector)) { # includes sliding of Semilandmarks
        if (is.null(outlines))
            stop("please specify outlines")
        dataslide <- Semislide(A, SMvector=SMvector,outlines=outlines,tol=tol,deselect=deselect,recursive=recursive,iterations=iterations,pairedLM=pairedLM,initproc=initproc, bending=bending,stepsize=stepsize)
        A <- dataslide
        for (i in 1:n)
            CS <- apply(A,3,cSize)
        if (CSinit==TRUE) { 
            for (i in 1:n)
                A[,,i] <- A[,,i]/CS[i]
        }
    }
###### create mirrored configs ######
    if (!is.null(pairedLM)) {
        Amir <- A
        for (i in 1:n) {
            Amir[,,i] <- A[,,i]%*%Mir
            Amir[c(pairedLM),,i] <- Amir[c(pairedLM[,2:1]),,i]
        }
        Aall <- bindArr(A,Amir,along=3)
    } else
        Aall <- A
    
###### proc fit of all configs ######
    cat("performing Procrustes Fit ")
    
    if (!is.null(use.lm)) { ### only use subset for rotation and scale
        proc <- ProcGPA(Aall[use.lm,,],scale=scale,CSinit=CSinit,reflection=reflect,pcAlign=pcAlign)
        tmp <- Aall
        for (i in 1:dim(Aall)[3]) {
            tmp[,,i] <- rotonmat(Aall[,,i],Aall[use.lm,,i],proc$rotated[,,i],scale=TRUE, reflection=reflect)
            if (center.part)
                tmp[,,i] <- scale(tmp[,,i], scale=FALSE) ## center shapes
            else
                orp <- FALSE
        }
        proc$rotated <- tmp
        proc$mshape <- arrMean3(tmp) ##calc new meanshape
    } else
        proc <- ProcGPA(Aall,scale=scale,CSinit=CSinit, reflection=reflect,pcAlign=pcAlign)
    
    procrot <- proc$rotated
    dimna <- dimnames(dataarray)
    if (!is.null(pairedLM))
        dimna[[3]] <- c(dimna[[3]],dimna[[3]])
    dimnames(proc$rotated) <- dimna
    meanshape <- proc$mshape
    rho <- NULL
    
    for (i in 1:n)
        rho[i] <- distfun(proc$rotated[,,i],proc$mshape)
    rmsrho <- sqrt(mean(rho^2))
    names(rho) <- dimnames(dataarray)[[3]]
    orpdata <- 0

###### project into tangent space ######
###test###        
    
    if (orp==TRUE && CSinit==TRUE)
        procrot <- orp(proc$rotated, mshape=proc$mshape)
    
    orpdata <- procrot
    dimnames(orpdata) <- dimna
    
###### calculate Symmetric means ######
    if (!is.null(pairedLM)) {
        ## generate symmetrized mean for each individual between original and mirrored configuration ###      		
        Symarray <- A
        for (i in 1:n)
            Symarray[,,i] <- (procrot[,,i]+procrot[,,n+i])/2
        ## generate deviation between each individual and its specific symmetrized mean ###    		
        Asymm <- A 
        for (i in 1:n)
            Asymm[,,i] <- (procrot[,,i]-Symarray[,,i])
        dimnames(Asymm) <- dimnames(dataarray)
    } else 
        Symarray <- procrot

    Symtan <- Symarray
    Symtan <- sweep(Symtan, 1:2, meanshape)
    tan <- vecx(Symtan)
    if (sizeshape) { 
        CSlog <- log(CS)-mean(log(CS))
        tan <- cbind(CSlog,tan)
    }
    dimnames(Symarray) <- dimnames(dataarray)
    
###### PCA Sym Component ###### 
    princ <- try(prcomp(tan),silent=TRUE)
    if (class(princ) == "try-error")
        princ <- eigenPCA(tan)

    values <- 0
    eigv <- princ$sdev^2
    values <- eigv[which(eigv > 1e-14)]
    lv <- length(values)
    PCs <- princ$rotation[,1:lv]
    PCscore_sym <- as.matrix(princ$x[,1:lv])
    rownames(PCscore_sym) <- dimnames(dataarray)[[3]]
    rownames(tan) <- rownames(PCscore_sym)

###### create a neat variance table for Sym ###### 
    if (length(values)==1){
        SymVar <- values
    } else {
        SymVar <- matrix(NA,length(values),3)
        SymVar[,1] <- values
        
        for (i in 1:length(values))
            SymVar[i,2] <- (values[i]/sum(values))*100
        SymVar[1,3] <- SymVar[1,2]
        for (i in 2:length(values))
            SymVar[i,3] <- SymVar[i,2]+ SymVar[i-1,3]
        colnames(SymVar) <- c("eigenvalues","% Variance","Cumulative %")
    }
    
###### PCA Asym Component ###### 
    asvalues <- 0
    PCs_Asym <- 0
    if (!is.null(pairedLM)) {
        asymmean <- arrMean3(Asymm)
        asymtan <- vecx(sweep(Asymm, 1:2, asymmean))[1:n,]
        pcasym <- try(prcomp(asymtan),silent=TRUE)
        if (class(pcasym) == "try-error")
            pcasym <- eigenPCA(asymtan)
        
        eigva <- pcasym$sdev^2
        asvalues <- eigva[which(eigva > 1e-14)]
        lva <- length(asvalues)
        PCs_Asym <- pcasym$rotation[,1:lva]
        PCscore_asym <- as.matrix(pcasym$x[,1:lva])
        rownames(PCscore_asym) <- dimnames(dataarray)[[3]]
        rownames(asymtan) <- rownames(PCscore_sym)
        
###### create a neat variance table for Asym ######
        if (length(asvalues)==1) {
            AsymVar <- asvalues
        } else {
            AsymVar <- matrix(NA,length(asvalues),3)
            AsymVar[,1] <- asvalues
            
            for (i in 1:length(asvalues))
                AsymVar[i,2] <- (asvalues[i]/sum(asvalues))*100
            
            AsymVar[1,3] <- AsymVar[1,2]
            for (i in 2:length(asvalues))
                AsymVar[i,3] <- AsymVar[i,2]+ AsymVar[i-1,3]
            
            colnames(AsymVar) <- c("eigenvalues","% Variance","Cumulative %")
        }
    }
    
###### output ######
    t1 <- Sys.time()
    cat(paste("Operation completed in",t1-t0,"secs\n"))
    if (!is.null(pairedLM)) {
        out <- (list(
            size=CS,rotated=proc$rotated[,,1:n],
            rotmir=proc$rotated[,,(n+1):(2*n)],
            Sym=Symarray,Asym=Asymm,asymmean=asymmean,
            mshape=(meanshape+asymmean), symmean=meanshape,
            Symtan=tan,Asymtan=asymtan,PCsym=PCs,PCscore_sym=PCscore_sym,
            eigensym=values,SymVar=SymVar,PCasym=PCs_Asym,
            PCscore_asym=PCscore_asym,eigenasym=asvalues,
            AsymVar=AsymVar,orpdata=orpdata[,,1:n],
            orpmir=orpdata[,,(n+1):(2*n)],
            rmsrho=rmsrho,rho=rho,dataslide= dataslide,
            pairedLM=pairedLM
            ))
        class(out) <- "symproc"
    } else {
        out <- (list(
            size=CS,rotated=proc$rotated,mshape=meanshape,tan=tan,PCs=PCs,
            PCscores=PCscore_sym,eigenvalues=values,Variance=SymVar,
            orpdata=orpdata[,,1:n] ,rmsrho=proc$rmsrho,rho=rho,
            dataslide= dataslide
            ))
        
        class(out) <- "nosymproc"
    }
    attributes(out) <- append(attributes(out),list(CSinit=CSinit,scale=scale,orp=orp,reflect=reflect))
    return(out)
    
}
#' @export       
print.nosymproc <- function(x,...) {
    cat(paste0(" No. of Specimens: ",dim(x$rotated)[3],"\n"))
        cat(paste0(" ",dim(x$rotated)[1]," Landmarks in ", dim(x$rotated)[2]," dimensions\n"))
    cat("\n Variance Table\n")
    print(as.data.frame(x$Var),row.names=FALSE)
}
    
#' @export       
print.symproc <- function(x,...) {
    cat(paste0(" No. of Specimens: ",dim(x$rotated)[3],"\n"))
    cat(paste0(" ",dim(x$rotated)[1]," Landmarks in ", dim(x$rotated)[2]," dimensions\n"))
    cat(paste0("    - of which there are ",dim(x$pairedLM)[1]," sets of paired Landmarks\n"))
    cat("\n Variance Table of Symmetric Component\n")
    print(as.data.frame(x$SymVar),row.names=FALSE)
    cat("\n Variance Table of Asymmetric Component\n")
    print(as.data.frame(x$AsymVar),row.names=FALSE)
}
    
