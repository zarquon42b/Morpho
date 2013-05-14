covDist <- function(s1,s2)
  {
    dims1 <- dim(s1);dims2 <- dim(s2)
    
    if (dims1[1] != dims1[2] || dims2[1] != dims2[2] || dims1[1] != dims2[1])
      stop("please provide covariance matrices with idential dimensionality")

    cdist <- sqrt(sum(log(eigen(solve(s1,s2))$values)^2))
    return(cdist)

  }
covPCA <- function(data,groups)
  {
    if (! is.factor(groups))
      groups <- as.factor(groups)
    lev <- levels(groups)
    nlev <- length(lev)
    covlist <- list()
    for (i in 1:nlev)
      covlist[[i]] <- cov(data[groups==lev[i],])

    V <- diag(0,nlev,nlev)
    for (i in 1:(nlev-1))
      {
        for (j in (i+1):(nlev))
          V[j,i] <- covDist(covlist[[j]],covlist[[i]])^2
      }
    V <- V+t(V)

    H <- matrix(-1/nlev,nlev,nlev)
    H <- H+diag(nlev)
    D <- (-1/2)*(H%*%V%*%H)
    eigenD <- eigen(D,symmetric = TRUE)
    eigenD$values <- eigenD$values[1:(nlev-1)]
    eigenD$vectors <- eigenD$vectors[,1:(nlev-1)]
    PCscores <- t(t(eigenD$vectors)*sqrt(eigenD$values))

    out <- list()
    out$PCscores <- PCscores
    out$Var <- eigenD$values/sum(eigenD$values)
    out$dist <- as.dist(V)
    out$eigen <- eigenD
    return(out)
  }

    
    
