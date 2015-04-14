createVarTable <- function (sdev, square = TRUE, rownames=NULL) {
    if (square) 
        sdev <- sdev^2
    sdsum <- sum(sdev)
    sdVar <- sdev/sdsum
    sdCum <- cumsum(sdVar)
    Variance <- data.frame(eigenvalues = sdev, exVar = sdVar, cumVar = sdCum)
    rownames(Variance) <- rownames
    return(Variance)
}
