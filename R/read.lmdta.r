#' read dta files
#' 
#' reads .dta files created by the software Landmark
#' http://graphics.idav.ucdavis.edu/research/EvoMorph
#' 
#' 
#' @param file a dta file
#' @param na specifies a value that indicates missing values
#' @return
#' \item{arr }{array containing landmarks dimnames will be Information of
#' ID and landmark names specified in Landmark}
#' \item{info }{Information extracted from the header of the dta file}
#' \item{idnames }{character vector containing the names of the individuals
#' as specified in the dta file}
#' 
#' @export
read.lmdta <- function(file="x", na=9999) {
    x <- file
    A <- readLines(x)
    em <- which(A=="")
    infoline <- grep("DIM=",A,ignore.case = TRUE)
    idlines <- infoline+min(which(A[-c(1:infoline)] != ""))
    endid <- idlines+min(which(A[-c(1:idlines)] == ""))
    idnames <- A[c((idlines):(endid-1))]
    info <- strsplit(A[infoline]," ")[[1]]
    n2 <- nchar(info[2])-1
    nspeci <- as.numeric(substr(info[2],1L,n2))
    ndim <- as.numeric(substr(info[6],5,nchar(info[6])))
    nlms <- as.numeric(info[3])/ndim
    eot <- endid
    B <- as.matrix(read.table(x,skip=eot),na.strings=as.numeric(info[5]))
    if (nrow(B) != nlms*nspeci)
        stop("number of landmarks in dataset not matching file header information")
    tt <- array(t(B),dim=c(ndim,nlms,nspeci))
    arr <- array(NA,dim=c(nlms,ndim,nspeci))
    for (i in 1:nspeci)
        arr[,,i] <- t(tt[,,i])

    nas <- which(arr == na)
    if (length(nas) > 0)
        arr[nas] <- NA
    
    dimnames(arr)[[3]] <- idnames
    return(list(arr=arr,info=info,idnames=idnames))
}
