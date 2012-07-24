read.mpp <- function(file)
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
    info <- lapply(data,subfuninfo)
    info <- gsub("\"","",unlist(lapply(info,function(x){x <- x[c(2,4)]})))
    info <- gsub("/>","",info)
    info <- matrix(info,length(points),2,byrow=T)
    infout <- data.frame(active=as.integer(info[,1]),name=info[,2])
    data <- lapply(data,subfun)
    tmp <- as.numeric(gsub("\"","",unlist(lapply(data,function(x){x <- x[c(2,4,6)]}))))
    tmp <- matrix(tmp,length(points),3,byrow=T)
    return(list(data=tmp,info=infout))
  }
