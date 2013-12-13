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
#' @export
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
#' @export
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
        if (!add)
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


