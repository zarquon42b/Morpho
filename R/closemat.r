closemeshKD <- function(x,mesh,k=50,sign=FALSE,cores=1,method=0,...)
    {
        if(.Platform$OS.type == "windows")
            cores <- 1
        if (is.null(mesh$normals))
            mesh <- adnormals(mesh)
        if (is.matrix(x)) {
            matr <- x
            x <- list()
            x$vb <- rbind(t(matr),1)
        } else {
            matr <- t(x$vb[1:3,])
        }
        vb <- (mesh$vb[1:3,])
        it <- (mesh$it)
        nvb <- ncol(vb)
        nit <- ncol(it)        
        nmat <- nrow(matr)
        dist <- rep(0,nmat)
        fptr <- dist
        bary <- barycenter(mesh)
        clostInd <- mcNNindex(bary,matr,k=k,cores=cores,...)
        storage.mode(k) <- "integer"
        clost <- c(0,0,0)
        storage.mode(fptr) <- "integer"
        storage.mode(it) <- "integer"
        storage.mode(nvb) <- "integer"    
        storage.mode(nit) <- "integer"
        storage.mode(method) <- "integer"
        storage.mode(nmat) <- "integer"
        storage.mode(matr) <- "double"
        storage.mode(vb) <- "double"
        storage.mode(dist) <- "double"
        storage.mode(clost) <- "double"
        storage.mode(clostInd) <- "integer"
        outmatr <- matr
        region <- fptr
        out <- .Fortran("matr_meshKD",ioMat=matr,nmat,vb,nvb,it,nit,clostInd,k,dist=dist,ptr=fptr,region,normals=mesh$normals[1:3,],sign,outnormals=t(matr),method)
        gc()
        x$vb[1:3,] <- t(out$ioMat)
        x$quality <- out$dist
        x$normals <- rbind(out$outnormals,1)
        x$ptr <- out$ptr
        class(x) <- "mesh3d"
        return(x)
    }
