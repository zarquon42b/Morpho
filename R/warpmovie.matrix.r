#' Creates a sequence of images showing predefined steps of warping two meshes
#' or landmark configurations (2D and 3D) into each other
#' 
#' given two landmark configurations or two meshes with the same amount of
#' vertices and faces (e.g a mesh and its warped counterpart), the starting
#' configuration/mesh will be subsequently transformed into the final
#' configuration/mesh by splitting the differences into a predefined set of
#' steps.
#' 
#' A series of png files will be saved to disk. These can be joined to animated
#' gifs by external programs such as imagemagick or used to create animations
#' in PDFs in a latex environment (e.g. latex package: aninmate).
#' @title Creates a sequence of images showing predefined steps of warping two meshes or landmark configurations (2D and 3D) into each other
#' 
#' @param x mesh to start with (object of class mesh3d)
#' @param y resulting mesh (object of class mesh3d), having the same amount of
#' vertices and faces than the starting mesh
#' @param n integer: amount of intermediate steps.
#' @param col color of the mesh
#' @param palindrome logical: if TRUE, the procedure will go forth and back.
#' @param folder character: output folder for created images (optional)
#' @param movie character: name of the output files
#' @param add logical: if TRUE, the movie will be added to the focussed
#' rgl-windows.
#' @param close logical: if TRUE, the rgl window will be closed when finished.
#' width and 200 the height of the image.
#' @param countbegin integer: number to start image sequence.
#' 
#' @param ask logical: if TRUE, the viewpoint can be selected manually.
#' @param radius numeric: define size of spheres (overides atuomatic size
#' estimation).
#' @param xland optional argument: add landmarks on mesh x
#' @param yland optional argument: add landmarks on mesh y
#' @param lmcol optional argument: color of landmarks xland and yland
#' @param links vector or list of vectors containing wireframe information to
#' connect landmarks (optional).
#' @param lwd numeric: controls width of lines defined by "links".
#' @param imagedim character of pattern "100x200" where 100 determines the
#' width and 200 the height of the image.
#' @param par list of graphial parameters: details can be found here:
#' \code{\link{par}}.
#' @param \dots additional arguments passed to \code{\link{shade3d}} (3D) or
#' \code{\link{points}} (2D).
#' @author Stefan Schlager
#' @seealso
#' \code{\link{ply2mesh},\link{file2mesh},\link{mesh2ply},\link{tps3d}}
#' 
#' @examples
#' 
#' 
#' ###3D example
#'  data(nose)##load data
#' \dontrun{
#' ##warp a mesh onto another landmark configuration:
#' warpnose.long <- tps3d(shortnose.mesh,shortnose.lm,longnose.lm)
#' 
#' warpmovie3d(shortnose.mesh,warpnose.long,n=15)## create 15 images.
#' 
#' ### ad some landmarks
#' warpmovie3d(shortnose.mesh,warpnose.long,n=15,xland=shortnose.lm,
#'             yland=longnose.lm)## create 15 images.
#'
#' 
#' ### restrict to landmarks
#' warpmovie3d(shortnose.lm,longnose.lm,n=15,movie="matrixmovie")## create 15 images.
#' 
#' ### the images are now stored in your current working directory and can
#' ### be concatenated to a gif using an external program such as
#' ### imagemagick.
#' }
#' ### 2D example
#' library(shapes)
#' bb <- procSym(gorf.dat)
#' ### morph superimposed first specimen onto sample mean
#' warpmovie2d(bb$rotated[,,1],bb$mshape,n=20,links=c(1,5,4:2,8:6,1),imagedim="600x400")
#' 
#' @rdname warpmovie3d
#' @export
warpmovie3d <- function (x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",...) UseMethod("warpmovie3d")

#' @rdname warpmovie3d
#' @method warpmovie3d matrix
#' @export
warpmovie3d.matrix <- function(x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",add=FALSE,close=TRUE,countbegin=0,ask=TRUE,radius=NULL,links=NULL,lwd=1,...)
{	#wdold <- getwd()
    if(!is.null(folder)) {
        if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") {
            folder <- paste(folder,"/",sep="")
            dir.create(folder,showWarnings=F)
            movie <- paste(folder, movie, sep="")
        }
    }
    if (!add)
        open3d()
    
    k <- dim(x)[1]
    if (is.null(radius))
        radius <- (cSize(x)/sqrt(k))*(1/80)
    ## get bbox
    bbox <- apply(rbind(x,y),2,range)
    bbox <- expand.grid(bbox[,1],bbox[,2],bbox[,3])
    points3d(bbox,col="white")
    
    for (i in 0:n) {
        mesh <- x
        mesh <- (i/n)*y+(1-(i/n))*x
        
        a <- spheres3d(mesh,col=col,radius=radius,...)
        if (!is.null(links)) {
            a1 <- lineplot(mesh,links,col=col,lwd=lwd)
            a <- append(a,a1)
        }
        if (i ==0 && ask==TRUE)
            readline("please select view and press return\n")
        
        filename <- sprintf("%s%04d.png", movie, countbegin+i)
        rgl.snapshot(filename,fmt="png")
        rgl.pop("shapes",id=a)
    }
    
    if (palindrome) {## go the other way ##
        for (i in 1:(n-1)) {
            mesh <- x
            mesh <- (i/n)*x+(1-(i/n))*y
            a <- spheres3d(mesh,col=col,radius=radius,...)
            if (!is.null(links)) {
                a1 <- lineplot(mesh,links,col=col,lwd=lwd)
                a <- append(a,a1)
            }
            filename <- sprintf("%s%04d.png", movie, countbegin+i+n)
            rgl.snapshot(filename,fmt="png")
            rgl.pop("shapes",id=a)	
        }
    }
    if (close)
        rgl.close()
}
#' @rdname warpmovie3d
#' @export
warpmovie2d <- function(x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",links=NULL,lwd=1,imagedim = "800x800",par=list(xaxt="n",yaxt="n",bty="n"),...)
{
    wdold <- getwd()
    widxheight <- as.integer(strsplit(imagedim, split = "x")[[1]])
    if(!is.null(folder)) {
        if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") {
            folder <- paste(folder,"/",sep="")
            dir.create(folder,showWarnings=F)
            setwd(folder)
        }
    }
    ##k <- dim(x)[1]
    ## get bbox
    bbox <- apply(rbind(x,y),2,range)
    bbox <- expand.grid(bbox[,1],bbox[,2])
    ## plot(bbox,col="white")
    
    for (i in 0:n) {
        mesh <- x
        mesh <- (i/n)*y+(1-(i/n))*x
        filename <- sprintf("%s%04d.png", movie, i)
        png(filename,width = widxheight[1], height = widxheight[2])
        par(par)
        plot(bbox,asp=1,cex=0,xlab="",ylab="")
        points(mesh,col=col,...)
        if (!is.null(links))
            lineplot(mesh,links,col=col,lwd=lwd)
        
        dev.off()
    }
    
    if (palindrome) {## go the other way ##
        for (i in 1:(n-1)) {
            mesh <- x
            mesh <- (i/n)*x+(1-(i/n))*y
            filename <- sprintf("%s%03d.png", movie, i+n)
            png(filename,width = widxheight[1], height = widxheight[2])
            par(par)
            plot(bbox,asp=1,cex=0)
            points(mesh,col=col,...)
            
            if (!is.null(links))
                lineplot(mesh,links,col=col,lwd=lwd)
            dev.off()
        }
    }
    setwd(wdold)
}
