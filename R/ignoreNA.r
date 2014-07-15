# find NAs to ignore missing data and recompute curves and outlines
ignoreNA <- function(datamatrix,outlines=NULL,SMvector=NULL,surp=NULL,pairedLM=0) {

    k <- nrow(datamatrix)
    ignore <- which(apply(datamatrix,1,function(x) x <- as.logical(sum(is.na(x)))))
    li <- length(ignore)
    lm.old <- c(1:k)[-ignore]
    mat.ptr <- matrix(c(1:(k-li),lm.old),k-li,2)
    ptr <- function(xo)	### define pointer function for indexing
        {
            if (length(which(ignore %in% xo))!= 0)
                xo <- xo[-which(xo %in% ignore)]
            for (i in 1:(k-li))
                xo[which(xo==mat.ptr[i,2])] <- mat.ptr[i,1]
            return(xo)
        }
    
    if (!is.null(outlines)) ### update outline indices
        outlines <- lapply(outlines,ptr)
    if (!is.null(surp)) 	### update surface indices
        surp <- ptr(surp)
    
    if (!is.null(SMvector)) ### of fixed/sliding definition
        SMvector <- ptr(SMvector)
    
    if (pairedLM[1]!=0){	### update paired landmarks indices
        count <- 0
        del <- NULL
        for (i in 1:dim(pairedLM)[1]) {	
            if (length(which(ignore %in% pairedLM[i,]))!=0) {
                count <- count+1
                del[count] <- i
            }
        }
        pairedLM <- pairedLM[-del,]
        if (is.vector(pairedLM))
            pairedLM <- t(as.matrix(pairedLM))
        
        if (dim(pairedLM)[1]==0) {
            pairedLM <- 0
        } else {
            pairedLM <- apply(pairedLM,2,ptr)
            if (is.vector(pairedLM))
                pairedLM <- t(as.matrix(pairedLM))
        }
    }
    datamatrix <- datamatrix[-ignore,]
    out <- list(datamatrix=datamatrix,ignore=ignore,outlines=outlines,SMvector=SMvector,surp=surp,pairedLM=pairedLM)
    return(out)
    
}


