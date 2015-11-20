#' Project semi-landmarks from a predefined atlas onto all specimen in a sample
#' 
#' Project semi-landmarks from a predefined atlas onto all specimen in a
#' sample. Various mechanisms are implemented to avoid errorneous placement on
#' the wrong surface layer (e.g. inside the bone).
#' 
#' This function allows the (relatively) easy projection of surface points
#' defined on an atlas onto all surface of a given sample by Thin-Plate Spline
#' deformation and additional mechanisms to avoid distortions. The algorithm
#' can be outlined as followed.  \enumerate{
#' \item relax curves (if specified) against atlas.
#' \item deform atlas onto targets by TPS based on predefined landmarks (and curves).
#' \item project coordinates on deformed atlas onto target mesh
#' \item 'inflate' or 'deflate' configuration along their normals to make sure
#' all coordinates are on the outside/inside
#' \item Project inflated points back onto surface along these normals.
#' \item Check if normals are roughly pointing into the same direction as those
#' on the (deformed) atlas.
#' \item Relax all points against atlas.
#' \item the predefined coordinates will note change afterwards!
#' 
#' }
#' 
#' @param atlas object of class "atlas" created by \code{\link{createAtlas}}
#' @param dat.array k x 3 x n array containing reference landmarks of the
#' sample or a matrix in case of only one target specimen.
#' @param path character: specify the directory where the surface meshes of the
#' sample are stored.
#' @param prefix character: prefix to the specimens names (stored in
#' \code{dimnames(dat.array)[[3]]}) to match the corresponding file names.
#' If \code{dat.array} has no dimnames (e.g. because it is a matrix - see example below),
#' this can also be a character vector
#' containing the filenames to which \code{fileext} will be appended.
#' @param fileext character: file extension of the surface meshes.
#' @param ray logical: projection will be along surface normals instead of
#' simple closest point search.
#' @param inflate inflate (or deflate - if negative sign) the semilandmarks
#' along the normals of the deformed atlas to make sure that they stay on the
#' outside (inside) of the target mesh.
#' @param tol numeric: threshold to follow the ray back after inflation. See
#' details below. If no surface is hit after \code{tol} mm, the simple closest
#' point will be used.
#' @param relax.patch logical: request relaxation minimising bending energy
#' toward the atlas.
#' @param keep.fix integer: rowindices of those landmarks that are not allowed
#' to be relaxed in case \code{relax.patch=TRUE}. If not specified, all
#' landmarks will be kept fix. This is preferably set during atlas creation with \code{createAtlas}:
#' In case you specified corrCurves on the atlas, you should define explicitly which landmarks
#' (also on the curves) are supposed to fix to prevent them from sliding.
#' @param rhotol numeric: maximum amount of deviation a hit point's normal is
#' allowed to deviate from the normal defined on the atlas. If
#' \code{relax.patch=TRUE}, those points exceeding this value will be relaxed
#' freely (i.e. not restricted to tangent plane).
#' @param silent logical: suppress messages.
#' @param mc.cores run in parallel (experimental stuff now even available on Windows).
#' On windows this will only lead to a significant speed boost for many configurations,
#' as all required packages (Morpho and Rvcg) need to be loaded by each newly spawned process.
#' @return array containing the projected coordinates appended to the
#' data.array specified in the input. In case dat.array is a matrix only a
#' matrix is returned.
#' @author Stefan Schlager
#' @seealso \code{\link{createAtlas}, \link{relaxLM}, \link{checkLM},
#' \link{slider3d}, \link{tps3d}}
#' @encoding utf8
#' @references Schlager S. 2013. Soft-tissue reconstruction of the human nose:
#' population differences and sexual dimorphism. PhD thesis,
#' \enc{Universitätsbibliothek}{Universitaetsbibliothek} Freiburg.  URL:
#' \url{http://www.freidok.uni-freiburg.de/volltexte/9181/}.
#' 
#' @examples
#' 
#' \dontrun{
#' data(nose)
#' require(rgl)
#' ###create mesh for longnose
#' longnose.mesh <- tps3d(shortnose.mesh,shortnose.lm,longnose.lm)
#' ## create atlas
#' fix <- c(1:5,20:21)
#' atlas <- createAtlas(shortnose.mesh, landmarks =
#'            shortnose.lm[fix,], patch=shortnose.lm[-c(1:5,20:21),])
#' ## view atlas
#' 
#' plotAtlas(atlas)
#' 
#' ## create landmark array with only fix landmarks
#' data <- bindArr(shortnose.lm[fix,], longnose.lm[fix,], along=3)
#' dimnames(data)[[3]] <- c("shortnose", "longnose")
#' 
#' ### write meshes to disk
#' mesh2ply(shortnose.mesh, filename="shortnose")
#' mesh2ply(longnose.mesh, filename="longnose")
#' 
#' patched <- placePatch(atlas, data, path="./", inflate=5)
#' ## now browse through placed patches
#' checkLM(patched, path="./", atlas=atlas)
#' 
#' ## same example with only one target specimen
#' data <- longnose.lm[fix, ]
#' 
#' patched <- placePatch(atlas, data, prefix="longnose", path="./", inflate=5)
#' wire3d(longnose.mesh,col=3)
#' spheres3d(patched)
#' }
#' @importFrom parallel makeCluster stopCluster
#' @export
placePatch <- function(atlas, dat.array, path, prefix=NULL, fileext=".ply", ray=TRUE, inflate=NULL,tol=inflate, relax.patch=TRUE, keep.fix=NULL, rhotol=NULL, silent=FALSE,mc.cores=1)
    {
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        if (!inherits(dat.array, "array") && !inherits(dat.array,"matrix"))
            stop("dat.array must be a numeric array or a matrix")
        if (is.null(keep.fix)) {
            if (is.null(atlas$keep.fix))
                keep.fix <- 1:dim(atlas$landmarks)[1]
            else
                keep.fix <- atlas$keep.fix
        }
        if (is.null(tol) && !is.null(inflate))
            tol <- inflate
        if (mc.cores > 1)
            silent <- TRUE
        
        patched <- place.patch(dat.array, path, atlas.mesh =atlas$mesh, atlas.lm = atlas$landmarks, patch =atlas$patch, curves=atlas$patchCurves, prefix=prefix, tol=tol, ray=ray, outlines=atlas$corrCurves, inflate=inflate, relax.patch=relax.patch, rhotol=rhotol, fileext=fileext, SMvector = keep.fix, silent=silent,mc.cores=mc.cores)
        return(patched)
    }
#' @importFrom foreach registerDoSEQ
place.patch <- function(dat.array,path,atlas.mesh,atlas.lm,patch,curves=NULL,prefix=NULL,tol=5,ray=T,outlines=NULL,SMvector=NULL,inflate=NULL,relax.patch=TRUE,rhotol=NULL,fileext=".ply", silent=FALSE, mc.cores=1)
    {
        if (.Platform$OS.type == "windows" && mc.cores > 1) {
            cl <- makeCluster(mc.cores)            
            registerDoParallel(cl=cl)
        } else if (mc.cores > 1) {
            registerDoParallel(cores = mc.cores)
        } else
            registerDoSEQ()
        
        k <- dim(dat.array)[1]
        deselect=TRUE
        fix <- which(c(1:k) %in% SMvector)
        patch.dim <- dim(patch)[1]
        usematrix <- FALSE
        if (! is.matrix(dat.array)) {
            n <- dim(dat.array)[3]
            out <- array(NA,dim=c((patch.dim+k),3,n))
            dimnames(out)[[3]] <- dimnames(dat.array)[[3]]
            name <- dimnames(dat.array)[[3]]
        } else {
            usematrix <- TRUE
            n <- 1
            out <- matrix(NA,(patch.dim+k),3)
            name <- NULL
        }
        nall <- patch.dim+k

        L <- CreateL(atlas.lm,output="Lsubk3")
        L1 <- CreateL(rbind(atlas.lm,patch),output="Lsubk3")
        meshpath <- paste(path,"/",prefix,name,fileext,sep="")
        i <- 0
        parfun <- function(i){
           
            tmp.name <- meshpath[i]
            tmp.mesh <- vcgImport(tmp.name)
            if (!usematrix)
                tmp.data <- projRead(dat.array[,,i],tmp.mesh,readnormals=TRUE)
            else
                tmp.data <- projRead(dat.array,tmp.mesh,readnormals=TRUE)
### relax existing curves against atlas ###
            if (!is.null(outlines)) {
                notinout <- which(! (1:k) %in% unlist(outlines))
                if (length(notinout))
                    SMvector <- unique(c(SMvector,notinout))
                sm <- SMvector
                deselcurve <- TRUE
                if (prod(length(unique(SMvector)) == k)) {
                    message("There are corresponding curves but no fix landmarks specified")
                    SMvector <- c(1:k)[which(! (1:k %in% unlist(outlines)))]
                    if (!length(SMvector)) {
                        SMvector <- 1:k
                        deselcurve <- FALSE
                    }
                }
                U <- .calcTang_U_s(t(tmp.data$vb[1:3,]),t(tmp.data$normals[1:3,]),SMvector=SMvector,outlines=outlines,surface=NULL,deselect=deselcurve)
                slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)
                tmp.data <- projRead(slide,tmp.mesh,readnormals=TRUE)
                tps.lm <- tps3d(patch,atlas.lm,slide)
            } else if (!is.null(SMvector) && is.null(outlines)) {
                sm <- SMvector
                tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            } else {
                sm <- 1:k
                tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            }
            
            slide <- t(tmp.data$vb[1:3,])
            slidenormals <- t(tmp.data$normals[1:3,])
            if (!usematrix)   #replace projected points with original for fix landmarks
                slide[fix,] <- dat.array[fix,,i]
            else
                slide[fix,] <- dat.array[fix,]

            if (!is.null(inflate) || !is.null(rhotol)) {
                atlas.warp <- tps3d(atlas.mesh,atlas.lm,slide, silent=silent)
                tps.lm <- projRead(tps.lm,atlas.warp,readnormals=TRUE,smooth=TRUE)
                warp.norm <- tps.lm$normals[1:3,]### keep projected normals
            }
### use for mullitlayer meshes to avoid projection inside
            if (!is.null(inflate)) {
                
                
                tps.lm$vb[1:3,] <- tps.lm$vb[1:3,]+inflate*tps.lm$normals[1:3,] ###inflate outward along normals
                tps.lm <- ray2mesh(tps.lm,tmp.mesh,inbound=TRUE,tol=tol) ### deflate in opposite direction
            } else {## just project warped patch on surface (suitable for singlelayer meshes)
                tps.lm <- projRead(tps.lm,tmp.mesh,readnormals=TRUE)
            }
            
            relax <- rbind(slide,t(tps.lm$vb[1:3,]))
            normals <- rbind(slidenormals,t(tps.lm$normals[1:3,]))
            surface <- c((k+1):(patch.dim+k))  ## define surface as appended to preset landmarks
            if (!is.null(curves)) {
                if (!is.list(curves))
                    curves <- list(curves)
                curves <- lapply(curves,function(x) x+k)
            }
            free <- NULL
### compare normals of projection and original points
            if (!is.null(rhotol)) {
                if (!relax.patch) {
                    relax.patch <- TRUE
                    cat("relaxing enabled because rhotol is set")
                }
                
                rho <- NULL
                rho <- sapply(1:patch.dim, function(j) {
                    out <- angle.calc(tps.lm$normals[1:3,j],warp.norm[1:3,j])
                    return(out)
                })
                rhoex <- which(rho > rhotol) 
                if (length(rhoex) > 0) {
                    free <- surface[rhoex]
                    surface <- surface[-rhoex]
                }
            }
            gc()
### end compare normals #### 
            
### relax patch against reference ###
            if (relax.patch){ ### relax against reference
                if (!is.list(outlines) && !is.null(outlines))
                    outlines <- list(outlines)
                outltmp <- append(outlines,curves) ## add curves from patch to predefined curves
                remout <- which(surface %in% unlist(curves))
               
                if (length(remout))
                    surface <- surface[-remout] ### remove patch curves from surface 
                if (!length(surface))
                    surface <- NULL
                
                U1 <- .calcTang_U_s(relax, normals,SMvector=sm,outlines=outltmp,surface=surface,free=free,deselect=deselect)
                tps.lm <- calcGamma(U1$Gamma0,L1$Lsubk3,U1$U,dims=3)[c((k+1):(patch.dim+k)),]
                tps.lm <- projRead(tps.lm,tmp.mesh,readnormals=FALSE)
            } else {# end relaxation ########################
                tps.lm <- t(tps.lm$vb[1:3,])
            }
            
            if (!usematrix)
                out <- rbind(dat.array[,,i],tps.lm)
            else
                out <- rbind(dat.array,tps.lm)
            return(out)
        }

        out <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.export=c("calcGamma",".calcTang_U_s"),.packages=c("Morpho","Rvcg")) %dopar% parfun(i)

        
        if (!usematrix && n > 1) {
            tmpout <- array(NA, dim=c(nall,3,n))
            for (i in 1:n) {
                if (is.matrix(out[[i]])) {
                    tmpout[,,i] <- out[[i]]
                } else {
                    warning(paste("matching for specimen",i,"failed with:",out[[i]]))
                }
            }
            out <- tmpout
            dimnames(out)[[3]] <-  dimnames(dat.array)[[3]]
        } else {
            out <- out[[1]]
            if (!is.matrix(out))
                stop("matching failed")
        }
        if (.Platform$OS.type == "windows" && mc.cores > 1)
            stopCluster(cl)
        return(out)
    }

