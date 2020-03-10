#' read fiducials from slicer4
#'
#' read fiducials from slicer4
#' @param x filename
#' @param na value to be replaced by NA
#' @return a k x 3 matrix with landmarks
#' @export
read.fcsv <- function(x,na=NULL) {
    raw <- readLines(x)
    getids <- grep("# columns",raw)
    myids <- raw[getids]
    myids <- gsub("# columns = ","",myids)
    myids <- unlist(strsplit(myids,split=","))

    xpos <- which(myids=="x")
    ypos <- which(myids=="y")
    zpos <- which(myids=="z")
    labelpos <- which(myids=="label")
    points <- which(!grepl("^#",raw))
    
    data <- strsplit(raw[points], split = ",")
    subfun <- function(x) {
        tmp <- strsplit(x[c(xpos,ypos,zpos)], split = "=")
        tmp <- unlist(tmp, recursive = F)
        return(tmp)
    }
    getnames <- function(x) {
        tmp <- strsplit(x[labelpos], split = "=")
        tmp <- unlist(tmp, recursive = F)
        return(tmp)
    }
    mynames <- unlist(lapply(data,getnames))
    data <- lapply(data, subfun)
    
    tmp <- as.numeric(unlist(data))
    tmp <- matrix(tmp, length(points), 3, byrow = T)
    if (!is.null(na)) {
       nas <- which(tmp == na)
       if (length(nas) > 0) 
           tmp[nas] <- NA 
    }
    rownames(tmp) <- mynames
    return(tmp)
}

#' write fiducials in slicer4 format
#'
#' write fiducials in slicer4 format
#' @param x matrix with row containing 2D or 3D coordinates
#' @param filename will be substituted with ".fcsv"
#' @param description optional: character vector containing a description for each landmark
#' @examples
#' require(Rvcg)
#' data(dummyhead)
#' write.fcsv(dummyhead.lm)
#' ## remove file
#' unlink("dummyhead.lm.fcsv")
#' @export 
write.fcsv <- function(x,filename=dataname,description=NULL) {
    dataname <- deparse(substitute(x))
    if (!grepl("*.fcsv$",filename))
        filename <- paste0(filename,".fcsv")
    cat("# Markups fiducial file version = 4.4\n# CoordinateSystem = 0\n# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n",file=filename)
    ptdim <- ncol(x)
    ptn <- nrow(x)
    if (ptdim == 2)
        x <- cbind(x,0)
    else if (ptdim != 3)
        stop("only 2D or 3D point clouds")
    rowids <- paste0("vtkMRMLMarkupsFiducialNode_",1:ptn)
    associatedNodeID <- rep("vtkMRMLVectorVolumeNode12",ptn)
    buffer <- rep(0,ptn)
    buffer <- cbind(buffer,0,0,1,1,1,0)
    desc <- rep("",ptn)
    names <- rownames(x)
    if (is.null(names))
        names <- paste0("F_",1:ptn)
    outframe <- data.frame(rowids,x,buffer,names,desc,associatedNodeID)
    write.table(format(outframe, scientific = F, 
            trim = T), file = filename, sep = ",", append = TRUE, 
            quote = FALSE, row.names = FALSE, col.names = FALSE, 
            na = "")
}
