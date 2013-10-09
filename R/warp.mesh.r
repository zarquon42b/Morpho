warp.mesh <- function(mesh,matr,matt,lambda=0,updateNormals=TRUE, silent=FALSE)
{
    vert <- t(mesh$vb[1:3,])
    if (!silent)
        cat("calculating spline...\n")
    warp <- tps3d(vert,matr,matt,lambda=lambda)
    mesh$vb <- rbind(t(warp),1)
    mesh$normals <- NULL
    testref <- rotonto(matr,matt)$reflect
    if(testref == 1)
        mesh <- conv2backf(mesh)
    
    if(updateNormals) {
        if (!silent)
            cat("updating normals...\n")
        mesh <- adnormals(mesh)
    }
    return(mesh)
}

