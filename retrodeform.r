#' @importFrom Rvcg vcgKDtree
wNeighbour <- function(mat,pairedLM,hmult=5) {
    npair <- nrow(pairedLM)
    P <- mat[pairedLM[,1],]
    Q <- mat[pairedLM[,2],]
    nnpd <- vcgKDtree(mat[pairedLM[,1],],mat[pairedLM[,1],],2)$distance[,-1]
   
    nnqd <- vcgKDtree(mat[pairedLM[,2],],mat[pairedLM[,2],],2)$distance[,-1]
    h <- hmult*mean(c(nnpd,nnqd))
    h2 <- h^2
    dp <- exp(-as.matrix(dist(mat[pairedLM[,1]])^2)/h2)
    dq <- exp(-as.matrix(dist(mat[pairedLM[,2]]))^2/h2)
    
    arr <- bindArr(dp,dq,along=3)
    PhiIJ <- apply(arr,1:2,min)
    diag(PhiIJ) <- 1
    getSi <- function(i) {
        Ptmp <- P*PhiIJ[,i]
        Qtmp <- Q*PhiIJ[,i]
        #Ptmp <- P
        #Qtmp <- Q
        trans <- colMeans(rbind(Qtmp,Ptmp))
        Psweep <- sweep(P,2,trans)
        Qsweep <- sweep(Q,2,trans)
        C <- 0
        for (j in 1:nrow(Ptmp))
            C <- PhiIJ[j,i]*(C+tcrossprod(Psweep[j,])+tcrossprod(Qsweep[j,]))
        C <- C/(2*npair-1)
        svdC <- eigen(C)
        Tinv <- svdC$vectors%*%diag(sqrt(svdC$values))%*%t(svdC$vectors)
        T <- svdC$vectors%*%diag(1/sqrt(svdC$values))%*%t(svdC$vectors)
        TP <- Psweep%*%T
        TQ <-Qsweep%*%T
        
        CQ <- 0
         for (j in 1:nrow(TP))
             CQ <- CQ+tcrossprod(TP[j,],TQ[j,])+ tcrossprod(TQ[j,],TP[j,])

        eigCQ <- eigen(CQ)
        Hstar <- Tinv%*%eigCQ$vectors
        Hstarunit <- apply(Hstar,2,function(x) x <- x/sqrt(sum(x^2)))
        w1 <- Hstar[,3]
        w1 <- w1/sqrt(sum(w1^2))
        n <- crossp(Hstar[,1],Hstar[,2])
        n <- n/sqrt(sum(n^2))
        if (crossprod(n,w1) < 0)
            n <- -n
        wtan <- tanplan(n)
        wtan <- cbind(wtan$z,wtan$y)
        m <- as.vector(wtan%*%t(wtan)%*%w1)
        m <- m/sqrt(sum(m^2))
        beta <- angle.calc(w1,-m)
        ny <- (w1-m)/2
        ny <- ny/sqrt(sum(ny^2))
        gamma <- (tan(beta/2))
        Si <- (gamma-1)*tcrossprod(ny)+diag(3)#
        ni <- as.vector(Si%*%w1)
        ni <- ni/sqrt(sum(ni^2))
        chk <- crossprod(ni,c(-1,0,0))
        if (chk < 0) {
            ni <- -ni
        }
        Qi <- t(rotonto(matrix(c(-1,0,0),1,3),matrix(ni,1,3),reflection = FALSE)$gam)
######
       
        return (list(Si=Si,ni=ni,Qi=Qi,PhiI=PhiIJ[,i]))
    }
    
    return(append(lapply(1:npair,getSi),list(allphi = PhiIJ)))
}

retrodeform <- function(mat,pairedLM,hmult=5,alpha=0.01) {
    #alpha <- 0.01
    npair <- nrow(pairedLM)
    P <- mat[pairedLM[,1],]
    Q <- mat[pairedLM[,2],]
    Pdiff <- lapply(1:npair,function(x){ out <- t(t(-P)+P[x,]);return(out )})
        Qdiff <- lapply(1:npair,function(x){ out <- t(t(-Q)+Q[x,]);return(out )})
    PQdiff <- lapply(1:npair,function(x){ out <- t(t(Q)-P[x,]);return(out )})
    QPdiff <- lapply(1:npair,function(x){ out <- t(t(-P)+Q[x,]);return(out )})
    precode <- wNeighbour(mat,pairedLM,hmult=hmult)
    sqrtPhiIJ <- (precode$allphi)
    diag(sqrtPhiIJ) <- 0
    Amat <- -8*(1+alpha)*sqrtPhiIJ
    diagA <- 8*(1+alpha)*colSums(sqrtPhiIJ)
    diag(Amat) <- diagA

    ## create Amat for case ri = si (x-dimension)
    ## first 2 terms
    Amatx <- -8*sqrtPhiIJ
    diag(Amatx) <- 8*colSums(sqrtPhiIJ)
    ## alpha terms
    Amatx1 <- precode$allphi*4*alpha
    newphi <- matrix(8,nrow(precode$allphi),ncol(precode$allphi))*alpha*precode$allphi
    diag(newphi) <- 16*alpha*diag(precode$allphi)
    diag(Amatx1) <- colSums(newphi)
    
    Amatx <- Amatx+Amatx1
    Bmat <- sqrtPhiIJ*2
    Bmat2 <- precode$allphi
    #Bmat <- Bmat*t(sqrt(precode$allphi))
    QiMiPx <- QiMiQx <-QiPx <- QiQx <- NULL
    QiMiPy <- QiMiQy <- QiPy <- QiQy <- NULL
    QiMiPz <- QiMiQz <- QiPz <- QiQz <- NULL
    for (i in 1:npair) {
        QiMiP <- t(precode[[i]]$Qi%*%precode[[i]]$Si%*%t(Pdiff[[i]]))
        QiMiQ <- t(precode[[i]]$Qi%*%precode[[i]]$Si%*%t(Qdiff[[i]]))
        QiP <- t(precode[[i]]$Qi%*%t(PQdiff[[i]]))
        QiQ <- t(precode[[i]]$Qi%*%t(QPdiff[[i]]))
        QiMiPx <- cbind(QiMiPx,(QiMiP[,1]))
        QiMiPy <- cbind(QiMiPy,(QiMiP[,2]))
        QiMiPz <- cbind(QiMiPz,(QiMiP[,3]))
        QiMiQx <- cbind(QiMiQx,(QiMiQ[,1]))
        QiMiQy <- cbind(QiMiQy,(QiMiQ[,2]))
        QiMiQz <- cbind(QiMiQz,(QiMiQ[,3]))
        QiPx <-  cbind(QiPx,(QiP[,1]))
        QiPy <-  cbind(QiPy,(QiP[,2]))
        QiPz <-  cbind(QiPz,(QiP[,3]))
        QiQx <-  cbind(QiQx,(QiQ[,1]))
        QiQy <-  cbind(QiQy,(QiQ[,2]))
        QiQz <-  cbind(QiQz,(QiQ[,3]))
    }
    bx <- 0
    by <- NULL
    bz <- NULL
    Bmatx <- precode$allphi*2
    diag(Bmatx) <- 4
    for (i in 1:npair) {
        tmp <- Bmat*QiMiPx
        tmpQ <- Bmat*QiMiQx
        tmpQi <- alpha*Bmatx*QiQx
        tmpPi <- alpha*Bmatx*QiPx
        bx[i] <- +sum(tmp[,i])-sum(tmp[i,])-sum(tmpQ[,i])+sum(tmpQ[i,])-sum(tmpPi[,i])-sum(tmpPi[i,])-sum(tmpQi[,i])-sum(tmpQi[i,])
       
        tmp <- Bmat*QiMiPy
        tmpQ <- Bmat*QiMiQy
        tmpQi <- alpha*Bmat*QiQy
        tmpPi <- alpha*Bmat*QiPy
        by[i] <- -sum(tmp[,i])+sum(tmp[i,])-sum(tmpQ[,i])+sum(tmpQ[i,])-sum(tmpPi[,i])+sum(tmpPi[i,])-sum(tmpQi[,i])+sum(tmpQi[i,])
        tmp <- Bmat*QiMiPz
        tmpQ <- Bmat*QiMiQz
        tmpQi <- alpha*Bmat*QiQz
        tmpPi <- alpha*Bmat*QiPz
        bz[i] <- -sum(tmp[,i])+sum(tmp[i,])-sum(tmpQ[,i])+sum(tmpQ[i,])-sum(tmpPi[,i])+sum(tmpPi[i,])-sum(tmpQi[,i])+sum(tmpQi[i,])
    }
    a <- cbind(Morpho:::armaGinv(Amatx)%*%(bx),Morpho:::armaGinv(Amat)%*%by,Morpho:::armaGinv(Amat)%*%bz)
    a1 <- a
    a1[,1] <- -a[,1]
    a <- rbind(a1,a)
    check <- as.logical(rotonto(mat,a)$reflect)
    if (check)
        a[,1] <- -a[,1]
    return(a)
}
