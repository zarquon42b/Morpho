#' import landmark data from csv files
#' 
#' import landmark data from csv files
#' 
#' 
#' @param file character: path to file containing landmark data.
#' @param x either a vector specifiing which rows are to be imported, or
#' character vector containing variable names to be sought for.
#' @param y a vector specifiing, which columns of the speradsheet ist to be
#' imported.
#' @param rownames integer: specifies columns, where variable names are stored.
#' @param header logical : if spreadsheet contains header-row.
#' @param dec character: defines decimal sepearator.
#' @param sep character: defines column seperator.
#' @return
#' \item{LM }{matrix containing imported data}
#' \item{NAs }{vector containing rows containing NAs}
#' @author Stefan Schlager
#' @seealso \code{\link{read.table}}
#' @export
readLandmarks.csv <- function(file, x, y=2:4, rownames=NULL, header=TRUE, dec=".", sep=";")
{	
	
    xlen <- length(x)
    ylen <- length(y)
    
    arr <- matrix(NA, xlen, ylen)
	if (is.factor(x))
            x <- as.character(x)
    
        if (is.character(x)) {### check if selection contains variable names 
            data <- read.table(file, header=header,dec=dec,sep=sep)
            dat <- NULL
            count <- 1
            if (is.null(rownames))
                stop("please specify column containing Landmark names!")
            
            rn <- data[,rownames]
            for (j in 1:length(x)) {
                check <- which(rn==x[j])
                if (length(check) == 0) {
                    warning(paste("dataset misses entry for Landmark",j))
                    data[99999, y] <- rep(NA,ylen)
                    dat[count] <- 99999
                }
                if (length(check) > 1) {
                    warning(paste("dataset contains landmark #",x[j],"with the same name - first match was used."))
                    dat[count] <- check[1]
                } else {
                    empty <- which(rn==x[j])
                    if (length(empty) !=0)
                        dat[count] <- which(rn==x[j])
                }
                count <- count+1
            }
            arr <- as.matrix(data[dat,y])
            rown <- x
        } else {
            data <- read.table(file,header=header,dec=dec,sep=sep)
            arr <- as.matrix(data[x,y])
            if (is.null(rownames))
                rown <- c(1:xlen)
            else
                rown <- data[x,rownames]
        }
    
    nas0 <- which(is.na(arr))	### check for NAs and store information about missing Landmark and individual
    nas1 <- as.integer(nas0/(xlen*ylen))+1
    nas <- nas1[-(which(duplicated(nas1)))]
    
    if (length(nas) > 0) {
        for (i in 1:length(nas)) {
            nas2 <- nas0[which(nas1==nas[i])]%%(xlen*ylen)
            nas2 <- nas2%%xlen
            nas2 <- nas2[-which(duplicated(nas2))]
            if (0 %in% nas2)
                {nas2[which(nas2==0)] <- xlen}
            if (length(nas2) > 0)
                nas2 <- sort(nas2)
        }
    } else 
        nas2 <- NULL
    
    dimnames(arr) <- list(rown,c("X","Y","Z"))
    return(list(LM=arr,NAs=nas2))
}
