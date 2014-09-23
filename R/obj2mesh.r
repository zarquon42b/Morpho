#' @rdname ply2mesh
#' @export
obj2mesh <- function(filename,adnormals=TRUE)
{	
    obj <- read.obj(filename)
    vert <- obj[which(obj[,1]=="v"),1:4]
    
    face <- obj[which(obj[,1]=="f"),1:4]
    vn <- obj[which(obj[,1]=="vn"),1:4]
    face.mat <- as.matrix(face[,2:4])
    vert.mat <- apply(vert[,2:4],2,as.numeric)
    
    if (length(grep("//",face.mat[1,1]))!=0) {	
        write.table(face.mat,file="facedump", quote = F, row.names = FALSE, col.names = FALSE, na = "",sep="//")
        face.mat <- read.table("facedump",sep="/")[,c(1,5,9)]	
        unlink("facedump")	
    }
    mesh <- list()
    class(mesh) <- "mesh3d"
    mesh$vb <- rbind(t(vert.mat),1)
    mesh$it <- t(face.mat)
    if (dim(vn)[1] != 0) {
        normals <- apply(vn[,2:4], 2, as.numeric)
        mesh$normals <- rbind(t(normals),1)
    }
    
    if (adnormals && is.null(mesh$normals))
        mesh <- vcgUpdateNormals(mesh)
    
    return(mesh)
}
