rot.proc <- function(x, y, scale=TRUE, weights=NULL, centerweight=FALSE, reflection=TRUE) {

    if (centerweight && !is.null(weights)) {
        xcent <- apply(x*weights,2,sum)
        ycent <- apply(y*weights,2,sum)
        x <- scale(x,scale=F,center=xcent)
        y <- scale(y,scale=F,center=ycent)
    }
### rotates 2 already centred matrices onto each other#
    if (!is.null(weights)) {
        Dn <- diag(weights)
        X1 <- Dn%*%x
        Y1 <- Dn%*%y
        XY <- crossprod(X1,Y1)
    } else
        XY <- crossprod(x,y)	

    sv1 <- svd(XY)
    gamm <- tcrossprod(sv1$v,sv1$u)
    if (sign(det(gamm)) < 1 && !reflection) {
        if (!reflection) {
            u <- sv1$u
            v <- sv1$v
            chk1 <- Re(prod(eigen(v)$values))
            chk2 <- Re(prod(eigen(u)$values))
            if ((chk1 < 0) && (chk2 > 0)) {
                v[, dim(v)[2]] <- v[, dim(v)[2]] * (-1)
                gamm <- v %*% t(u)
            }
            if ((chk2 < 0) && (chk1 > 0)) {
                u[, dim(u)[2]] <- u[, dim(u)[2]] * (-1)
                gamm <- v %*% t(u)
            }
        }
    }
    del <- sv1$d
    ctrace <- function(MAT) sum(diag(crossprod(MAT)))
    if (scale) {
       
        if (!is.null(weights))
            bet <- sum(del)/ctrace(Y1)
        else
            bet <- sum(del)/ctrace(y)
        yrot <- bet*y%*%gamm}
    else
        yrot <- y%*%gamm
    return(yrot)
}
