CreateL <- function(matrix,lambda=0)
{
    if (dim(matrix)[2] ==3) {
        k <- dim(matrix)[1]
        q1 <- matrix(c(rep(1,k)),k,1)
        K <- matrix(0,k,k)
        Q <- cbind(q1,matrix)
        O <- matrix(c(rep(0,16)),4,4)
        
        storage.mode(K) <- "double"
        storage.mode(matrix) <- "double"
        K <- .Fortran("createL",K,nrow(K),matrix,ncol(matrix))[[1]]
        
        diag(K) <- lambda
        L <- rbind(cbind(K,Q),cbind(t(Q),O))
        L1 <- try(solve(L),silent=TRUE)
        if (class(L1)=="try-error") {
            cat("CreateL: singular matrix: general inverse will be used.\n")
            L1 <- ginv(L)		
        }
        Lsubk <- L1[1:k,1:k]
        Lsubk3 <- matrix(0,3*k,3*k)
        Lsubk3[1:k,1:k] <- Lsubk
        Lsubk3[(k+1):(2*k),(k+1):(2*k)] <- Lsubk
        Lsubk3[(2*k+1):(3*k),(2*k+1):(3*k)] <- Lsubk;
        
        return(list(L=L,Linv=L1,Lsubk=Lsubk,Lsubk3=Lsubk3))
    }
    else if (dim(matrix)[2] == 2) {
        
        out <- CreateL2D(matrix, lambda)
        return(out)
    } else
        stop("only works for matrices with 2 or 3 columns")
}
CreateL2D <- function(matrix, lambda=0)
{
    k <- dim(matrix)[1]
    q1 <- matrix(c(rep(1,k)),k,1)
    K <- matrix(0,k,k)
    Q <- cbind(1,matrix)
    O <- matrix(0,3,3)

    for (i in 1:k) {
        for (j in 1:k) {
            r2 <- sum((matrix[i,]-matrix[j,])^2)
            K[i,j] <- r2*log(r2)
        }
    }
    K[which(is.na(K))] <- 0
    diag(K) <- lambda
    L <- rbind(cbind(K,Q),cbind(t(Q),O))
    
	L1 <- try(solve(L),silent=TRUE)
    	if (class(L1)=="try-error") {
            cat("singular matrix: general inverse will be used.\n")
            L1 <- ginv(L)		
        }
    Lsubk <- L1[1:k,1:k]
    Lsubk3 <- rbind(cbind(Lsubk,matrix(0,k,k)),cbind(matrix(0,k,k),Lsubk))
    return(list(L=L,Linv=L1,Lsubk=Lsubk,Lsubk3=Lsubk3))
}
