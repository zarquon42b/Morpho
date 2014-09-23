#' Import 3D surface mesh files
#' 
#' imports 3D mesh files and store them as an R .object of class mesh3d
#' 
#'  
#' @title Import 3D surface mesh files
#' @param filename character: path to file
#' @param adnormals Logical: If the file does not contain normal information,
#' they will be calculated in R: Can take some time.
#' @param readnormals Logical: Import vertex normals (if available), although
#' no face information is present.
#' @param readcol Logical: Import vertex colors (if available).
#' @param clean Logical: Delete dumpfiles.
#' @param silent logical: suppress messages.
#' @return
#' \item{mesh }{list of class mesh3d - see rgl manual for further details,
#' or a matrix containing vertex information or a list containing vertex and
#' normal information}
#' 
#' @examples
#' 
#' data(nose)
#' mesh2ply(shortnose.mesh)
#' mesh <- ply2mesh("shortnose.mesh.ply")
#' 
#' mesh2obj(shortnose.mesh)
#' mesh2 <- obj2mesh("shortnose.mesh.obj")
#' @rdname ply2mesh
#' @export
ply2mesh <- function (filename, adnormals = TRUE,readnormals=FALSE,readcol=FALSE, silent=FALSE)
{
    x <- filename
    A <- readLines(x, n = 100)
    end <- which(A == "end_header")
    infos <- A[1:end]
    vertinfo <- strsplit(A[grep("element vertex", infos)], " ")
    faceinfo <- strsplit(A[grep("element face", infos)], " ")
    colmat <- NULL
    material <- NULL
    ##if (length(grep("property list uchar float texcoord",A))==1) 
    qualinfo <- grep("property float quality", infos)
    vertbegin <- grep("element vertex",infos)
    facebegin <- grep("element face",infos)
    if (length(qualinfo)==1) {
        qualine <- qualinfo-vertbegin
        qual <- TRUE
    } else {
        qual <- FALSE
    }
    fn <- as.numeric(faceinfo[[1]][3])
    vn <- as.numeric(vertinfo[[1]][3])
    ##vert.all <- read.table(x, skip = end, sep = " ", nrows = vn,colClasses="numeric")
    vert.all <- scan(x, skip = end, nlines=vn,quiet=TRUE)
    vert.all <- matrix(vert.all,vn,length(vert.all)/vn,byrow=T)
    vert <- vert.all[, 1:3]
    vert.n <- NULL
    quality <- NULL

    if (length(grep("property float nx", infos)) == 1) {
        normstart <- grep("property float nx", infos)-vertbegin
        vert.n <- t(vert.all[, normstart:(normstart+2)])
        if (qual) {
            quality <- as.vector(vert.all[,qualine])
        }
    }

    if (readcol) {##check for colored vertices
        color <- grep("property uchar red",infos)
        if (length(color > 0)) {
            if (color[1] < facebegin) {
                colbegin <- color[1]-vertbegin
                colmat <- vert.all[,colbegin:(colbegin+2)]
                colmat <- rgb(colmat[,1],colmat[,2],colmat[,3],maxColorValue=255)
            }
        }
    }
#### read face data
    datatype <- integer()
    if (length(grep("float",A[facebegin:end])) > 0)
        datatype <- double()
    if (fn !=0) {
        face.all <- scan(x, skip = end+vn, nlines=fn,quiet=TRUE,what=datatype)
        face.all <- matrix(face.all,fn,length(face.all)/fn,byrow=T)
        ##face.all <- read.table(x, skip = end + vn, nrows = fn,colClasses="integer")
        face <- t(face.all[, 2:4]+1)
        
        if (!is.null(colmat)) {
            #dimface <- dim(face)
            colfun <- function(x)
                {
                    x <- colmat[x]
                    return(x)
                }
            material$color <- matrix(colfun(face),dim(face))
        }
        mesh <- list(vb = rbind(t(vert), 1), it = face, primitivetype = "triangle", material = material,normals = vert.n)      
        class(mesh) <- c("mesh3d", "shape3d")
    } else {
        if (is.null(vert.n) || readnormals==FALSE) {
            if (!silent)
                cat(paste("mesh contains no faces. Vertices will be stored in a",vn,"x 3 matrix\n"))
            mesh <- vert
        } else if (readnormals) {
            if (!silent)
            cat(paste("mesh contains no faces. vertices and vertex normals are stored in a list\n"))
            mesh <- list(vb=t(vert),normals=vert.n)
            class(mesh) <- "mesh3d"
        }	
    }
### generate object of class mesh3d ###	
    
### add TexCoords ###
    if (length(grep("property list uchar float texcoord",A))==1 && length(grep("comment TextureFile",A))==1) {
        texn <- face.all[1,5]
        tex <- face.all[,c(6:(6+(texn-1)))]
        mesh$tex <- t(tex)
        mesh$TextureFile <- strsplit(A[grep("comment TextureFile", infos)], " ")[[1]][3]
    } 
    
### check for normals and update if required ###
    if (fn !=0)	{
        if (adnormals && "mesh3d" %in% class(mesh)) {
            if (is.null(mesh$normals)) {
                if (!silent)
                    cat("calculating normals...\n")
                mesh <- vcgUpdateNormals(mesh)
            }
        }
    }
    if (qual && !is.matrix(mesh))
        mesh$quality <- quality

    return(mesh)
}

