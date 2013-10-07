warp.mesh <- function(mesh,matr,matt,lambda=0,updateNormals=TRUE)
{
    vert <- t(mesh$vb[1:3,])
    cat("calculating spline...\n")
    warp <- tps3d(vert,matr,matt,lambda=lambda)
    mesh$vb <- rbind(t(warp),1)
    mesh$normals <- NULL
    testref <- rotonto(matr,matt)$reflect
    if(testref == 1)
        mesh <- conv2backf(mesh)
    
    if(updateNormals) {
        cat("updating normals...\n")
        mesh <- adnormals(mesh)
    }
    return(mesh)
}

