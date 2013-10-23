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
