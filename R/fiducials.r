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
#' @param slicer4.11 logical: Slicer changed their fiducial format in version >= 4.11. Set TRUE if you use the latest Slicer version
#' @examples
#' require(Rvcg)
#' data(dummyhead)
#' write.fcsv(dummyhead.lm)
#' ## remove file
#' unlink("dummyhead.lm.fcsv")
#' @export 
write.fcsv <- function(x,filename=dataname,description=NULL,slicer4.11=FALSE) {
    dataname <- deparse(substitute(x))
    if (!grepl("*.fcsv$",filename))
        filename <- paste0(filename,".fcsv")
    if (!slicer4.11)
        cat("# Markups fiducial file version = 4.4\n# CoordinateSystem = 0\n# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n",file=filename)
    else
        cat("# Markups fiducial file version = 4.11\n# CoordinateSystem = LPS\n# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n",file=filename)
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

#' convert data from LPS to RAS space and back
#'
#' convert data from LPS to RAS space and back
#'
#' @param x mesh of class \code{mesh3d} or a matrix with 3D Landmarks
#'
#' @return returns the rotated data
#' @details As e.g. the Slicer versions >= 4.11 are using LPS space, it might be needed to transform data like fiducials and surface models from and back to that space.
#' @export
LPS2RAS <- function(x) {
    x <- applyTransform(x,diag(c(-1,-1,1,1)))
    return(x)
}

#' read Landmarks from Slicer in Json format
#'
#' read Landmarks from Slicer in Json format
#' @param x path to json file
#'
#' returns matrix or list of matrices with imported landmark coordinates
#' @importFrom jsonlite read_json
#' @export
read.slicerjson <- function(x) {
    
    
    mydata <- read_json(x,T)$markups
    message(paste0("Importing: ",mydata$type, "\nCoordinate System: ", mydata$coordinateSystem))
    
    
    cp <- mydata$controlPoints
    helpfun <- function(z) {
        mat <- t(sapply(z$position,rbind))
        labels <- z$label
        rownames(mat) <- labels
        return(mat)
    }
    
    cp <- lapply(cp,helpfun)
    cp <- lapply(cp,function(x){ attributes(x) <- append(attributes(x),list(coordinateSystem=mydata$coordinateSystem));return(x)})
    if (length(cp) == 1)
        cp <- cp[[1]]

   
    return(cp)
}
