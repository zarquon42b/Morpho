mirrMissingArray <- function(x,pairedLM) {
    n <- dim(x)[3]
    out <- x
    for (i in 1:n) {
        out[,,i] <- mirrMissing(x[,,i],pairedLM=pairedLM)
    }
    return(out)
}
        
mirrMissing <- function(x,pairedLM) {
    data <- x
    
    checkvec <-NULL
    count <- 0
    k <- nrow(x)
    checklist <- NA
    unilatNA <- NULL
    for (j in 1:k) {
        if (TRUE %in% is.na(data[j,])) {
            count <- count+1
            checklist[count] <- j
            checkvec <- 1
        }
    }
    affected <- affectCol <- goodPaired <- NULL
    for (i in 1:nrow(pairedLM)) {
        if (prod(is.na(x[pairedLM[i,],]))) {
            message("one landmark of each side must be present")
            unilatNA <- append(unilatNA,pairedLM[i,])
            
        }
        checkPaired <- is.na(x[pairedLM[i,],])
        if (TRUE %in% checkPaired) {
            affected <- append(affected,i)
            checkCol <- (apply(checkPaired,1,prod))
            affectCol <- append(affectCol,which(!as.logical(checkCol)))
            goodPaired <- append(goodPaired,pairedLM[i,which(!as.logical(checkCol))])
        }
    }
    if (is.null(affected)) {
        message("no missing bilateral landmarks")
        return(x)
    }
   
    if (!prod(checklist %in% pairedLM)){
        warning("missing landmarks are not bilateral")
        unilatNA <- append(unilatNA,checklist[which(! checklist %in% pairedLM)])
    }
     
    
    xmir <- x %*% diag(c(-1,1,1))#mirror landmarks
    xmir[c(pairedLM),] <- xmir[c(pairedLM[,2:1]),]##relabel landmarks
    xref <- xmir[-c(unilatNA,pairedLM[affected,]),]
    xtar <- x[-c(unilatNA,pairedLM[affected,]),]
    xout <- tps3d(xmir,xref,xtar)
    xout[goodPaired,] <- x[goodPaired,]
    return(xout)

}
