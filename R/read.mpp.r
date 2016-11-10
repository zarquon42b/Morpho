#' Read saved pick-points from meshlab
#' 
#' imports pick points selected with meshlab
#' 
#' 
#' @param file file to import
#' @param info logical: if TRUE, addtional infos are returned
#' @return if \code{info=FALSE}:
#' 
#' a matrix containing picked-points coordinates (only those tagged as active).
#' 
#' if \code{info=TRUE}: a list containing
#' \item{data }{matrix containing coordinates - including points tagged as inactive}
#' \item{info }{additional info contained in file.}
#' @author Stefan Schlager
#' @seealso \code{\link{read.pts}}
#' 
#' @export
read.mpp <- function(file, info=FALSE) {
    raw <- readLines(file)
    points <- grep("<point ",raw)
    data <- strsplit(raw[points],split="\ ")

    getcoord <- function(x) {
        xc <- grep("x=",x)
        yc <- grep("y=",x)
        zc <- grep("z=",x)
        namec <- grep("name=",x)
        activec <- grep("active=",x)
        dataind <- c(xc,yc,zc)
        infoind <- c(namec,activec)
        return(list(dataind=dataind,infoind=infoind))
    }
    
    getinds <- lapply(data,getcoord)
    datamat <- infomat <- NULL
    for (i in 1:length(getinds)) {
        tmppoints <- data[[i]][getinds[[i]]$dataind]
        tmppoints <- as.numeric(gsub("[a-z,\\,\",\\=,\\/,\\>]","",tmppoints))
        datamat <- rbind(datamat,tmppoints)
        tmpinfo <- data[[i]][getinds[[i]]$infoind]
        tmpinfo <- gsub("[a-z,\\,\",\\=,\\/,\\>]","",tmpinfo)
        infomat <- rbind(infomat,tmpinfo)
    }
    rownames(infomat) <- NULL
    colnames(infomat) <- c("name","active")
    infomat <- as.data.frame(infomat)
    rownames(datamat) <- infomat$name
    
    if (info) 
        return(list(data=datamat,info=infomat))
    else {
        datamat <- datamat[which(infomat$active == 1),]
        return(datamat)
    }
}


