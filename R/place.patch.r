#' Create an atlas needed in placePatch
#' 
#' Create an atlas needed in placePatch
#' 
#' 
#' @param mesh triangular mesh representing the atlas' surface
#' @param landmarks matrix containing landmarks defined on the atlas, as well
#' as on each specimen in the corresponding sample.
#' @param patch matrix containing semi-landmarks to be projected onto each
#' specimen in the corresponding sample.
#' @param corrCuves integer vector specifiyng the rowindices of
#' \code{landmarks} to be curves defined on the atlas AND each specimen.
#' @param patchCurves integer vector specifiyng the rowindices of \code{patch}
#' to be curves only defined on the atlas.
#' @return Returns a list of class "atlas".  Its content is corresponding to
#' argument names.
#' @note This is a helper function of \code{\link{placePatch}}.
#' @seealso \code{\link{placePatch}, \link{plotAtlas}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data(nose)
#' atlas <- createAtlas(shortnose.mesh, landmarks =
#'             shortnose.lm[c(1:5,20:21),], patch=shortnose.lm[-c(1:5,20:21),])
#' 
#' @export createAtlas
createAtlas <- function(mesh, landmarks, patch, corrCuves=NULL, patchCurves=NULL)
    {
        atlas <- list()
        class(atlas) <- "atlas"
        atlas$mesh <- mesh
        atlas$landmarks <- landmarks
        atlas$patch <- patch
        atlas$corrCurves <- corrCuves
        atlas$patchCurves<- patchCurves
        return(atlas)
    }


#' visualize an atlas defined by createAtlas
#' 
#' visualize an atlas defined by createAtlas
#' 
#' If \code{legend=TRUE}, a plot with a legend will open where coloring of the
#' 3D-spheres is specified.
#' 
#' @param atlas object of class atlas created by \code{\link{createAtlas}}.
#' @param pt.size size of plotted points/spheres. If \code{point="s"}.
#' \code{pt.size} defines the radius of the spheres. If \code{point="p"} it
#' sets the variable \code{size} used in \code{point3d}.
#' @param alpha value between 0 and 1. Sets transparency of mesh 1=opaque 0=
#' fully transparent.
#' @param render if \code{render="w"}, a wireframe will be drawn, if
#' \code{render="s"}, the mesh will be shaded.
#' @param point how to render landmarks. "s"=spheres, "p"=points.
#' @param meshcol color to render the atlas mesh
#' @param add logical: if TRUE, a new rgl window is opened.
#' @param legend logical: request plot of legend specifying landmark coloring.
#' @param cols vector containing colors for each coordinate type cols[1]=landmarks, cols[2]=patch, cols[3]=corrCurves, cols[4]=patchCurves.
#' @return returns invisible vector containing \code{rgl.id} of rendered
#' objects.
#' @seealso \code{\link{placePatch}, \link{createAtlas}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data(nose)
#' atlas <- createAtlas(shortnose.mesh, landmarks =
#'            shortnose.lm[c(1:5,20:21),], patch=shortnose.lm[-c(1:5,20:21),])
#' plotAtlas(atlas)
#' 
#' @export plotAtlas
plotAtlas <- function(atlas, pt.size=NULL, alpha=1, render=c("w","s"), point=c("s", "p"), meshcol="white", add=TRUE, legend=TRUE,cols=2:5)
    {
        outid <- NULL
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        point <- point[1]
        ## set point/sphere sizes
        radius <- pt.size
        if (is.null(radius)) {
            if (point == "s")
                radius <- (cSize(atlas$landmarks)/sqrt(nrow(atlas$landmarks)))*(1/30)
            else
                radius <- 10
        }
        size <- radius
        render <- render[1]
        if (point == "s") {
            rendpoint <- spheres3d
        } else if (point == "p") {
            rendpoint <- points3d
        } else {
            stop("argument \"point\" must be \"s\" for spheres or \"p\" for points")
        }
        if (render=="w") {
            rend <- wire3d
        } else {
            rend <- shade3d
        }
        if (add)
            open3d()
        if (!is.null(atlas$mesh))
            outid <- rend(atlas$mesh, col=meshcol, alpha=alpha)
        ## plot reference landmarks and patch
        landm <- atlas$landmarks
        if (!is.null(atlas$corrCurves))
            landm <- landm[-unlist(atlas$corrCurves),]
        patch <- atlas$patch
         if (!is.null(atlas$patchCurves))
            patch <- patch[-unlist(atlas$patchCurves),]
        outid <- c(outid, rendpoint(landm, col=cols[1], radius=radius, size=size))
        outid <- c(outid,rendpoint(patch,col=cols[2],radius=radius/2, size=size/2))
        ## plot reference curves
        if (!is.null(atlas$corrCurves))
            outid <- c(outid, rendpoint(atlas$landmarks[unlist(atlas$corrCurves),],col=cols[3],radius=radius, size=size))
                           
        if (!is.null(atlas$patchCurves))
            outid <- c(outid,rendpoint(atlas$patch[unlist(atlas$patchCurves),],col=cols[4],radius=radius/2,size=size/2))
        if (legend) {
            plot(0,0, xlab="", ylab="", axes =F, cex=0,xlim=c(-1,1), ylim=c(-1,1))
            legend(-1,1, pch=20, cex=2, col=cols, legend=c("landmarks", "patch", "curves on all specimen", "curves only on atlas"))
        }
        invisible(outid)
    }



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
#' \item deform atlas onto targets by TPS based on predefined landmarks (and curves).
#' \item project coordinates on deformed atlas onto target mesh
#' \item 'inflate' or 'deflate' configuration along their normals to make sure
#' all coordinates are on the outside/inside
#' \item Project inflated points back onto surface along these normals.
#' \item Check if normals are roughly pointing into the same direction as those
#' on the (deformed) atlas.
#' \item Relax all points against atlas.
#' 
#' }
#' 
#' @param atlas object of class "atlas" created by \code{\link{createAtlas}}
#' @param dat.array k x 3 x n array containing reference landmarks of the
#' sample or a matrix in case of only one target specimen.
#' @param path character: specify the directory where the surface meshes of the
#' sample are stored. In case \code{dat.array} is a matrix, \code{path} is
#' needed to specifiy not only the directory, but the specific surface mesh
#' file (without fileextension if specified by \code{fileext}). See examples
#' below.
#' @param prefix character: append a prefix to the specimens names (stored in
#' \code{dimnames(dat.array)[[3]]}) to match the corresponding file names.
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
#' landmarks will be kept fix.
#' @param rhotol numeric: maximum amount of deviation a hit point's normal is
#' allowed to deviate from the normal defined on the atlas. If
#' \code{relax.patch=TRUE}, those points exceeding this value will be relaxed
#' freely (i.e. not restricted to tangent plane).
#' @param silent logical: suppress messages.
#' @return array containing the projected coordinates appended to the
#' data.array specified in the input. In case dat.array is a matrix only a
#' matrix is returned.
#' @note needs additional command line tools "trimesh-tools" installed
#' (\url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/}).
#' @author Stefan Schlager
#' @seealso \code{\link{createAtlas}, \link{relaxLM}, \link{checkLM},
#' \link{slider3d}, \link{warp.mesh}}
#' @encoding utf8
#' @references Schlager S. 2013. Soft-tissue reconstruction of the human nose:
#' population differences and sexual dimorphism. PhD thesis,
#' \enc{Universitätsbibliothek}{Universitaetsbibliothek} Freiburg.  URL:
#' \url{http://www.freidok.uni-freiburg.de/volltexte/9181/}.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{
#' data(nose)
#' require(rgl)
#' ###create mesh for longnose
#' longnose.mesh <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
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
#' patched <- placePatch(atlas, data, path="longnose", inflate=5)
#' wire3d(longnose.mesh,col=3)
#' spheres3d(patched)
#' }
#' 
#' @export placePatch
placePatch <- function(atlas, dat.array, path, prefix=NULL, fileext=".ply", ray=TRUE, inflate=NULL,tol=inflate, relax.patch=TRUE, keep.fix=NULL, rhotol=NULL, silent=FALSE)
    {
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        if (is.null(keep.fix))
            keep.fix <- 1:dim(atlas$landmarks)[1]
        if (is.null(tol) && !is.null(inflate))
            tol <- inflate
        
        patched <- place.patch(dat.array, path, atlas.mesh =atlas$mesh, atlas.lm = atlas$landmarks, patch =atlas$patch, curves=atlas$patchCurves, prefix=prefix, tol=tol, ray=ray, outlines=atlas$corrCurves, inflate=inflate, relax.patch=relax.patch, rhotol=rhotol, fileext=fileext, SMvector = keep.fix, silent=silent)
        return(patched)
    }

place.patch <- function(dat.array,path,atlas.mesh,atlas.lm,patch,curves=NULL,prefix=NULL,tol=5,ray=T,outlines=NULL,SMvector=NULL,inflate=NULL,relax.patch=TRUE,rhotol=NULL,fileext=".ply", silent=FALSE)
    {
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
        
        L <- CreateL(atlas.lm)
        L1 <- CreateL(rbind(atlas.lm,patch))
        
        for(i in 1:n) {
            tmp.name <- paste(path,prefix,name[i],fileext,sep="")
            if (!usematrix)
                tmp.data <- projRead(dat.array[,,i],tmp.name,readnormals=TRUE, ignore.stdout=silent)
            else
                tmp.data <- projRead(dat.array,tmp.name,readnormals=TRUE, ignore.stdout=silent)
            
### relax existing curves against atlas ###
            if (!is.null(outlines)) {
                sm <- SMvector
                U <- .calcTang_U_s(t(tmp.data$vb[1:3,]),t(tmp.data$normals[1:3,]),SMvector=SMvector,outlines=outlines,surface=NULL,deselect=deselect)
                slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
                tmp.data <- projRead(slide,tmp.name,readnormals=TRUE, ignore.stdout=silent)
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
                                      
### use for mullitlayer meshes to avoid projection inside
            if (!is.null(inflate)) {
                atlas.warp <- warp.mesh(atlas.mesh,atlas.lm,slide, silent=silent)
                tps.lm <- projRead(tps.lm,atlas.warp,readnormals=TRUE,smooth=TRUE, ignore.stdout = silent)
                warp.norm <- tps.lm$normals[1:3,]### keep projected normals
                
                tps.lm$vb[1:3,] <- tps.lm$vb[1:3,]+inflate*tps.lm$normals[1:3,] ###inflate outward along normals
                tps.lm <- ray2mesh(tps.lm,tmp.name,inbound=TRUE,tol=tol,angmax=rhotol, ignore.stdout=silent) ### deflate in opposite direction
            } else {## just project warped patch on surface (suitable for singlelayer meshes)
                tps.lm <- projRead(tps.lm,tmp.name,readnormals=TRUE, ignore.stdout=silent )
            }
            
            relax <- rbind(slide,t(tps.lm$vb[1:3,]))
            normals <- rbind(slidenormals,t(tps.lm$normals[1:3,]))
            
            surface <- c((k+1):(patch.dim+k))  ## define surface as appended to preset landmarks
            free <- NULL
### compare normals of projection and original points
            if (!is.null(rhotol)) {
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
                outltmp <- append(outlines,curves) ## add curves from patch to predefined curves
                remout <- which(surface %in% curves)
                
                if (length(remout) > 0)
                    surface <- surface[-remout] ### remove patch curves from surface 
                if (length(surface)==0)
                    surface <- NULL
                
                U1 <- .calcTang_U_s(relax, normals,SMvector=sm,outlines=outltmp,surface=surface,free=free,deselect=deselect)
                tps.lm <- calcGamma(U1$Gamma0,L1$Lsubk3,U1$U,dims=3)$Gamatrix[c((k+1):(patch.dim+k)),]
                tps.lm <- projRead(tps.lm,tmp.name,readnormals=FALSE, ignore.stdout = silent)
            } else {# end relaxation ########################
                tps.lm <- t(tps.lm$vb[1:3,])
            }
            
            if (!usematrix)
                out[,,i] <- rbind(dat.array[,,i],tps.lm)
            else
                out <- rbind(dat.array,tps.lm)
        }
        return(out)
    }

