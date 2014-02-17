#' Read saved pick-points from meshlab
#' 
#' imports pick points selected with meshlab
#' 
#' 
#' @param file file to import
#' @param info logical: if TRUE, addtional infos are returned
#' @return if \code{info=FALSE}:
#' 
#' a matrix containing picked-points coordinates
#' 
#' if \code{info=TRUE}: a list containing
#' \item{data }{matrix containing coordinates}
#' \item{info }{additional info contained in file}
#' @author Stefan Schlager
#' @seealso \code{\link{read.pts}}
#' 
#' @export
read.mpp <- function(file, info=FALSE)
  {
    raw <- readLines(file)
    points <- grep("point x",raw)
    data <- strsplit(raw[points],split="\ ")

    subfun <- function(x)
      {
        tmp <- strsplit(x[3:5],split="=")
        tmp <- unlist(tmp,recursive=F)
        
        return(tmp)
      }
    subfuninfo <- function(x)
      {
        tmp <- strsplit(x[6:7],split="=")
        tmp <- unlist(tmp,recursive=F)
        
        return(tmp)
      }
    infoin <- lapply(data,subfuninfo)
    infoin <- gsub("\"","",unlist(lapply(infoin,function(x){x <- x[c(2,4)]})))
    infoin <- gsub("/>","",infoin)
    infoin <- matrix(infoin,length(points),2,byrow=T)
    infout <- data.frame(active=as.integer(infoin[,1]),name=infoin[,2])
    data <- lapply(data,subfun)
    tmp <- as.numeric(gsub("\"","",unlist(lapply(data,function(x){x <- x[c(2,4,6)]}))))
    tmp <- matrix(tmp,length(points),3,byrow=T)
    rownames(tmp) <- infout$name
    tmp <- tmp[which(infout$active == 1),]
    if (info)
        return(list(data=tmp,info=infout))
    else
        return(tmp)
  }
