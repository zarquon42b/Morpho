#' export mesh objects to disk.
#' 
#' export an object of class \code{mesh3d} or a set of coordinates to a common
#' mesh file.
#' 
#' 
#' @title export mesh objects to disk
#' @param x object of class \code{mesh3d} - see rgl documentation for further
#' details or a matrix containing vertices, this can either be a \code{k x 3}
#' or a \code{3 x k} matrix, with rows or columns containing vertex
#' coordinates.
#' @param filename character: Path/name of the requested output - extension
#' will be added atuomatically. If not specified, the file will be named as the
#' exported object.
#' @param col Writes color information to ply file. Can be either a single
#' color value or a vector containing a color value for each vertex of the
#' mesh.
#' @param writeNormals logical: if TRUE, existing normals of a mesh are written
#' to file - can slow things down for very large meshes.
#' @author Stefan Schlager
#' @note meshes containing quadrangular faces will be converted to triangular meshes by splitting the faces.
#' @seealso \code{\link{ply2mesh}, \link{quad2trimesh} }
#' 
#' @examples
#' 
#' require(rgl)
#' vb <- c(-1.8,-1.8,-1.8,1.0,1.8,-1.8,-1.8,1.0,-1.8,1.8,-1.8,1.0,1.8,
#' 1.8,-1.8,1.0,-1.8,-1.8,1.8,1.0,1.8,
#' -1.8,1.8,1.0,-1.8,1.8,1.8,1.0,1.8,1.8,1.8,1.0)
#' it <- c(2,1,3,3,4,2,3,1,5,5,7,3,5,1,2,2,6,5,6,8,7,7,5,6,7,8,4,4,3,7,4,8,6,6,2,4)
#' vb <- matrix(vb,4,8) ##create vertex matrix
#' it <- matrix(it,3,12) ## create face matrix
#' cube<-list(vb=vb,it=it)
#' class(cube) <- "mesh3d"
#'\dontrun{
#' shade3d(cube,col=3) ## view the green cube
#' }
#' mesh2ply(cube,filename="cube") # write cube to a file called cube.ply
#' @rdname mesh2ply
#' @export
mesh2ply <- function(x, filename=dataname, col=NULL, writeNormals=FALSE)
{	
    dataname <- deparse(substitute(x))
    if (is.matrix(x)) {
        dimsx <- dim(x)
        if (dimsx[2] == 3 && dimsx[1] != 3)
            x <- t(x)
        x <- list(vb=x)
    }
    if (!is.null(x$ib))
        x <- quad2trimesh(x)
    
    filename <- paste(filename,".ply",sep="")
    vert <- x$vb[1:3,]
    vert <- round(vert,digits=6)
    if (!is.null(x$it)) {
        face <- x$it-1
        fd <- 3
        fn <- dim(face)[2]
    } else
        fn <- 0
    
    vert.all <- vert
    vn <- dim(vert)[2]
    ##vn.all <- 3
    texfile <- x$TextureFile

    if (is.null(col) && !is.null(x$material$color)) {
        col=rep("#FFFFFF",vn)
        tmp1 <- data.frame(it=as.vector(x$it))
        tmp1$rgb <- as.vector(x$material$color)
        tmp1 <- unique(tmp1)
        col[tmp1$it] <- tmp1$rgb
    }
    if (!writeNormals)
        x$normals <- NULL
    
### start writing to file ###
    
### write header ###
    cat("ply\nformat ascii 1.0\ncomment MORPHO generated\n",file=filename)
    
### check for Texture information and write to header ###	
    
    if (is.character(texfile))
        cat(paste("comment TextureFile ",texfile,"\n",sep=""),file=filename,append=TRUE) 
    
    
### write vertex infos to header ###
    
    v.info <- paste("element vertex ",vn,"\n","property float x\nproperty float y\nproperty float z\n",sep="")
    
    cat(v.info,file=filename,append=TRUE)	
    if (!is.null(x$normals)) {
        cat("property float nx\nproperty float ny\nproperty float nz\n",file=filename,append=TRUE)
        norma <- round(x$normals[1:3,],digits=6)
        vert.all <- rbind(vert,norma)	
        ##vn.all <- 6		
    }
    if (!is.null(col))    
        v.info <- cat("property uchar red\nproperty uchar green\nproperty uchar blue\n",file=filename,append=T)
    
### write face infos and texture infos to header and finish ###	
    
    cat(paste("element face ",fn,"\n",sep=""),file=filename,append=TRUE)	
    if(!is.null(x$tex) && is.character(texfile)) {
        cat("property list uchar int vertex_indices\nproperty list uchar float texcoord\nend_header\n",file=filename,append=TRUE)	
    } else
        cat("property list uchar int vertex_indices\nend_header\n",file=filename,append=TRUE)	
    
### write vertices and vertex normals ###
    vert.all <- data.frame(t(vert.all))
    if (!is.null(col)) {
        if (is.vector(col)) {
            colout <- t(col2rgb(col))
        } else
            colout <- matrix(col2rgb(col),vn,3,byrow=T)
        
        vert.all <- cbind(vert.all,colout)
    }
    
    write.table(format(vert.all,scientific=F,trim=T),file=filename,sep=" ",append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")	
    ##write face and Texture information ###	
    if (!is.null(x$it)){
	if(!is.null(x$tex) && !is.null(x$TextureFile)) {
            tex <- t(x$tex)
            texn <- dim(tex)[2]
            faceraw <- rbind(fd,face)
            facetex <- t(cbind(texn,tex))
            write.table(format(t(rbind(faceraw,facetex)),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")			
        } else 
            write.table(format(t(rbind(fd,face)),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")	
    }
}
