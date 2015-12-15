#' calculates and visualises distances between surface meshes or 3D coordinates
#' and a surface mesh.
#' 
#' calculates the distances from a mesh or a set of 3D coordinates to another
#' at each vertex; either closest point or along the normals
#' 
#' this function needs the command line tools from the Auxiliaries section in
#' \url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries}
#' installed.
#' 
#' @title  calculates and visualises distances between surface meshes or 3D coordinates and a surface mesh.
#' @param x reference mesh; object of class "mesh3d" or a n x 3 matrix
#' containing 3D coordinates.
#' @param mesh2 target mesh: either object of class "mesh3d" or a character
#' pointing to a surface mesh (ply, obj or stl file)
#' @param distvec vector: optional, a vector containing distances for each
#' vertex of mesh1, if distvec != NULL, x will be ignored.
#' @param from numeric: minimum distance to be colorised; default is set to 0
#' mm
#' @param to numeric: maximum distance to be colorised; default is set to the
#' maximum distance
#' @param steps integer: determines break points for color ramp: n steps will
#' produce n-1 colors.
#' @param ceiling logical: if TRUE, the next larger integer of "to" is used
#' @param rampcolors character vector: specify the colors which are used to create a colorramp.
#' @param NAcol character: specify color for values outside the range defined by \code{from} and \code{to}.
#' @param file character: filename for mesh and image files produced. E.g.
#' "mydist" will produce the files mydist.ply and mydist.png
#' @param imagedim character of type 100x200 where 100 determines the width and
#' 200 the height of the image.
#' @param uprange numeric between 0 and 1: restricts "to" to a quantile of
#' "to", if to is NULL.
#' @param ray logical: if TRUE, the search is along vertex normals.
#' @param raytol maximum distance to follow a normal.
#' @param raystrict logical: if TRUE, only outward along normals will be sought for closest points.
#' @param save logical: save a colored mesh.
#' @param plot logical: visualise result as 3D-plot and distance charts
#' @param sign logical: request signed distances. Only meaningful, if mesh2 is
#' specified or distvec contains signed distances.
#' @param tol numeric: threshold to color distances within this threshold
#' green.
#' @param displace logical: if TRUE, displacement vectors between original and
#' closest points are drawn colored according to the distance.
#' @param shade logical: if FALSE, the rendering of the colored surface will be
#' supressed.
#' @param method accepts: "vcglib" and "morpho" (and any abbreviation). vcglib
#' uses a command line tool using vcglib headers, morpho uses fortran routines
#' based on a kd-tree search for closest triangles.
#' @param type character: "s" shows coordinates as spheres, while "p" shows 3D
#' dots.
#' @param radius determines size of spheres; if not specified, optimal radius
#' size will be estimated by centroid size of the configuration.
#' @param add logical: if TRUE, visualization will be added to the rgl window currently in focus
#' @param scaleramp logical: if TRUE, the colorramp will be symmetrical for signed distances: spanning from \code{-max(from,to)} to \code{max(from,to)}.
#' @param \dots additional arguments passed to \code{\link{shade3d}}. See
#' \code{\link{rgl.material}} for details.
#' @return Returns an object of class "meshDist" if the input is a surface mesh
#' and one of class "matrixDist" if input is a matrix containing 3D
#' coordinates.
#' \item{colMesh }{object of mesh3d with colors added}
#' \item{dists }{vector with distances}
#' \item{cols }{vector with color values}
#' \item{params }{list of parameters used}
#' @author Stefan Schlager
#' @seealso \code{\link{render.meshDist}}, , \code{\link{export.meshDist}},
#' \code{\link{shade3d}}
#' @references Detection of inside/outside uses the algorithm proposed in:
#' 
#' Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance
#' Fields From Triangle Meshes. Informatics and Mathematical Modelling, .
#' 
#' @examples
#' 
#' data(nose)##load data
#' ##warp a mesh onto another landmark configuration:
#' warpnose.long <- tps3d(shortnose.mesh, shortnose.lm, longnose.lm)
#' \dontrun{
#' mD <- meshDist(warpnose.long, shortnose.mesh)
#' ##now change the color ramp
#' render(mD,rampcolors = c("white","red"))
#' }
#' #use unsigned distances and a ramp from blue to red
#' #color distances < 0.01 green:
#' \dontrun{
#' meshDist(warpnose.long, shortnose.mesh, rampcolors = c("blue", "red"),sign=FALSE, tol=0.5)
#' }
#' @rdname meshDist
#' @export
meshDist <- function(x,...) UseMethod("meshDist")

#' @rdname meshDist
#' @method meshDist mesh3d
#' @importFrom Rvcg vcgClostKD
#' @importFrom colorRamps blue2green2red
#' @export
meshDist.mesh3d <- function(x, mesh2=NULL, distvec=NULL, from=NULL, to=NULL, steps=20, ceiling=FALSE,  rampcolors=colorRamps::blue2green2red(steps-1),NAcol="white", file="default", imagedim="100x800", uprange=1, ray=FALSE, raytol=50, raystrict=FALSE, save=FALSE, plot=TRUE, sign=TRUE, tol=NULL, displace=FALSE, shade=TRUE, method=c("vcglib", "morpho"), add=FALSE,scaleramp=TRUE,...)
  {
    method=substring(method[1],1L,1L)
    neg=FALSE
    NAcol <- colorRampPalette(NAcol)(1)
    #ramp <- blue2green2red(steps-1)
    ramp <- colorRampPalette(rampcolors)(steps-1)
    if (is.null(distvec)) {
        if(!ray) {
            if (method == "v") {
                promesh <- vcgClostKD(x,mesh2,sign=T)
            } else {
                promesh <- closemeshKD(x,mesh2,sign=T)
            }
            clost <- promesh$vb[1:3,]
            dists <- promesh$quality
            distsOrig <- dists
            if (!sign)
                dists <- abs(dists)
        } else {
            promesh <- ray2mesh(x,mesh2,tol=raytol,mindist=!raystrict)
            clost <- promesh$vb[1:3,]
            dists <- promesh$distance
            distsOrig <- dists
            if (!sign)
              dists <- abs(dists)
        }
    } else {
        clost <- NULL
        dists <- distvec
        distsOrig <- dists
        if (!sign)
            dists <- abs(dists)
    }  
    
    if (is.null(from)) {
        mindist <- min(dists)
        if (sign && mindist < 0 ) {
            from <- quantile(dists,probs=(1-uprange)) 
            neg <- TRUE            
        } else {
            from <- 0
        }
    }
    if (from < 0)
        neg <- TRUE
    if (is.null(to))
        to <- quantile(dists,probs=uprange)    
    if(ceiling)
        to <- ceiling(to)
    
    to <- to+1e-10
    colseq <- seq(from=from,to=to,length.out=steps)
    coldif <- colseq[2]-colseq[1]
    if (neg && sign) {
        negseq <- length(which(colseq<0))
        poseq <- steps-negseq
        maxseq <- max(c(negseq,poseq))
        ramp <- colorRampPalette(rampcolors)(maxseq*2)
        if (scaleramp) {
                ramp <- colorRampPalette(rampcolors)(maxseq*2)
                ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
            }
            else
                ramp <- colorRampPalette(rampcolors)(steps-1)
        distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
                                        #distqual[which(distqual < 1)] <- steps+10
    } else if (from > 0) {
          distqual <- ceiling(((dists-from)/coldif)+1e-14)
      } else {
            distqual <- ceiling((dists/coldif)+1e-14)
        }
    distqual[which(distqual < 1)] <- steps+10
    colorall <- ramp[distqual]
    
    if (!is.null(tol)) {
        if (sign) {
            tol <- c(-tol,tol)
        } else {
            tol <- c(0,tol)
        }
        good <- which(abs(dists) < tol[2])
        colorall[good] <- "#00FF00"
    }   
    
    colfun <- function(x){x <- colorall[x];return(x)}
    x$material$color <- matrix(colfun(x$it),dim(x$it))
    x$material$color[is.na(x$material$color)] <- NAcol
    colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
    params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol,rampcolors=rampcolors,NAcol=NAcol,scaleramp=scaleramp)
    out <- list(colMesh=x,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
    class(out) <- "meshDist"

    if (plot)
        render(out,output=FALSE,displace=displace,shade=shade,add=add, ...)
    if (save)
        export(out,file=file,imagedim=imagedim)
    invisible(out)
}



#' plot or save the results of meshDist
#' 
#' Visualise or save the results of meshDist to disk.
#' 
#' render.meshDist renders the colored mesh and displays the color ramp and
#' returns an object of class "meshDist".  export.meshDist exports the colored
#' mesh as ply file and the color chart as png file.
#' 
#' @title plot or save the results of meshDist
#' @param x object of class meshDist
#' @param from numeric: minimum distance to color; default is set to 0 mm
#' @param to numeric: maximum distance to color; default is set to the maximum
#' distance
#' @param steps integer: determines how many intermediate colors the color ramp
#' has.
#' @param ceiling logical: if TRUE, the next larger integer of "to" is used
#' @param uprange numeric between 0 and 1: restricts "to" to a quantile of
#' "to", if to is NULL.
#' @param tol numeric: threshold to color distances within this threshold
#' green.
#' @param rampcolors character vector: specify the colors which are used to create a colorramp.
#' @param NAcol character: specify color for values outside the range defined by \code{from} and \code{to}.
#' @param displace logical: if TRUE, displacement vectors between original and
#' closest points are drawn colored according to the distance.
#' @param shade logical: if FALSE, the rendering of the colored surface will be
#' supressed.
#' @param sign logical: request signed distances to be visualised.
#' @param file character: filename for mesh and image files produced. E.g.
#' "mydist" will produce the files mydist.ply and mydist.png
#' @param imagedim character of pattern "100x200" where 100 determines the
#' width and 200 the height of the image.
#' @param type character: "s" shows coordinates as spheres, while "p" shows 3D
#' dots.
#' @param radius determines size of spheres; if not specified, optimal radius
#' size will be estimated by centroid size of the configuration.
#' @param add logical: if TRUE, visualization will be added to the rgl window currently in focus
#' @param scaleramp if TRUE the ramp colors get scaled symmetrically into positive and negative direction.
#' @param \dots for render.meshDist: additional arguments passed to
#' \code{\link{shade3d}}. See \code{\link{rgl.material}} for details.
#' @author Stefan Schlager
#' @seealso \code{\link{meshDist}}, \code{\link{shade3d}}
#' 
#' @rdname render
#' @export
#'
render <- function(x,...) UseMethod("render")

#' @rdname render
#' @method render meshDist
#' @export
render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,shade=TRUE,sign=NULL,add=FALSE,scaleramp=NULL,...)
  {
    clost <- x$clost
    dists <- x$dists
    distsOrig <- dists
    colorall <- x$cols
    colramp <- x$colramp
    params <- x$params
    distqual <- x$distqual
    
    if (!add) {
        if (rgl.cur() !=0)
            rgl.clear()
    }
    if (!is.null(from) || !is.null(to) || !is.null(uprange) ||  !is.null(tol)  ||  !is.null(sign) || !is.null(steps) || !is.null(rampcolors) || !is.null(NAcol) || !is.null(scaleramp)) {
        neg=FALSE
        colMesh <- x$colMesh
        if(is.null(steps))
            steps <- x$params$steps
        if (is.null(rampcolors))
            rampcolors <- x$params$rampcolors
        if (is.null(NAcol))
            NAcol <- x$params$NAcol
        if(is.null(sign))
          sign <- x$params$sign
        if (!sign) {
            distsOrig <- dists
            dists <- abs(dists)
        }
        if(is.null(ceiling))
            ceiling <- x$params$ceiling
        if(is.null(uprange))
          uprange <- x$params$uprange
                 
        if (is.null(from)) {
            mindist <- min(dists)
            if (sign && mindist < 0 ) {
                from <- quantile(dists,probs=(1-uprange)) 
                neg <- TRUE            
            } else {
                from <- 0
            }             
        }
        if (is.null(scaleramp))
            scaleramp <- x$scaleramp
        if (from < 0)
            neg <- TRUE
        if (is.null(to))
            to <- quantile(dists,probs=uprange)    
        if(ceiling)
            to <- ceiling(to)
        
        to <- to+1e-10
        #ramp <- blue2green2red(maxseq*2)
        ramp <- colorRampPalette(rampcolors)(steps-1)
        colseq <- seq(from=from,to=to,length.out=steps)
        coldif <- colseq[2]-colseq[1]
        if (neg && sign) {
            
            negseq <- length(which(colseq<0))
            poseq <- steps-negseq
            maxseq <- max(c(negseq,poseq))
                                        
            if (scaleramp) {
                ramp <- colorRampPalette(rampcolors)(maxseq*2)
                ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
                
            }
            else
                ramp <- colorRampPalette(rampcolors)(steps-1)
            distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
            #distqual[which(distqual < 1)] <- steps+10
        } else if (from > 0) {
              distqual <- ceiling(((dists-from)/coldif)+1e-14)
          } else {
                distqual <- ceiling((dists/coldif)+1e-14)
            }
        distqual[which(distqual < 1)] <- steps+10
        colorall <- ramp[distqual]
        if (!is.null(tol)) {
            if (sign) {
                tol <- c(-tol,tol)
            } else {
                tol <- c(0,tol)
            }
            good <- which(abs(dists) < tol[2])
            colorall[good] <- "#00FF00"
        }
        colfun <- function(x){x <- colorall[x];return(x)}
        colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
        colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
        #colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
        colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
    } else {
        if (is.null(tol))
            tol <- x$params$tol
        
        colramp <- x$colramp
        colMesh <- x$colMesh
    }
    if (shade)
        shade3d(vcgUpdateNormals(colMesh),specular="black",...)
    if (displace) {
        dismesh <- colMesh
        vl <- dim(colMesh$vb)[2]
        dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
        dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
        dismesh$material$color <- rbind(colorall,colorall,colorall)
        wire3d(dismesh,lit=FALSE)
    }
    diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
    image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
    if (!is.null(tol)) {
        if (sum(abs(tol)) != 0)
            image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col="green",useRaster=TRUE,add=TRUE)
    }
    params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol,rampcolors=rampcolors,NAcol=NAcol)
    out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
                                
    class(out) <- "meshDist"
    invisible(out)
}
