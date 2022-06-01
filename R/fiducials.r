#' read fiducials from slicer4
#'
#' read fiducials from slicer4
#' @param x filename
#' @param na value to be replaced by NA
#' @param lps2ras logical: if the coordinate system is LPS and \code{lps2ras=TRUE}, the data will be rotated into the RAS space by inverting the first two dimensions using \code{\link{LPS2RAS}}.
#' @return a k x 3 matrix with landmarks
#' @export
read.fcsv <- function(x,na=NULL,lps2ras=FALSE) {
    raw <- readLines(x)
    getids <- grep("# columns",raw)
    myids <- raw[getids]
    myids <- gsub("# columns = ","",myids)
    myids <- unlist(strsplit(myids,split=","))
    LPS=FALSE
    getCooType <- grep("# CoordinateSystem",raw)
    coo <- gsub("# CoordinateSystem = ","",raw[getCooType])
    if (coo == "LPS")
        LPS=TRUE
    
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
    if (LPS && lps2ras)
        tmp <- LPS2RAS(tmp)
    if (!is.null(na)) {
        nas <- which(abs(tmp) == na)
        if (length(nas) > 0) 
            tmp[nas] <- NA 
    }
    
    if (nrow(tmp) == length(mynames))
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
#' @param lps2ras logical: if the coordinate system is LPS and \code{lps2ras=TRUE}, the data will be rotated into the RAS space by inverting the first two dimensions using \code{\link{LPS2RAS}}.
#' @param na value to be replaced by NA
#' @return returns matrix or list of matrices with imported landmark coordinates
#' @importFrom jsonlite read_json
#' @export
read.slicerjson <- function(x,lps2ras=FALSE,na=NULL) {
    
    
    mydata <- read_json(x,T)$markups
    message(paste0("Importing: ",mydata$type, "\nCoordinate System: ", mydata$coordinateSystem))
    LPS=FALSE
    if (mydata$coordinateSystem == "LPS")
        LPS = TRUE
    
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

     if (LPS && lps2ras)
         cp <- LPS2RAS(cp)

    if (!is.null(na)) {
        nas <- which(abs(cp) == na)
        if (length(nas) > 0) 
            cp[nas] <- NA 
    }
    return(cp)
}



#' Export landmarks (or any 3D coordinates) to the new slicer json format
#'
#' Export landmarks (or any 3D coordinates) to the new slicer json format
#' @param x k x 3 matrix containing 3D coordinates
#' @param filename will be substituted with ".mrk.json"
#' @param type character: specify type of coordinates. Can be any of "Fiducial", "Curve", "ClosedCurve".
#' @param coordinateSystem character: specify coordinate system the data are in. Can be "LPS" or "RAS".
#' @param labels character or character vector containing landmark labels. 
#' 
#' @importFrom jsonlite write_json
#' @export
write.slicerjson <- function(x,filename=dataname,type=c("Fiducial","Curve","ClosedCurve"),coordinateSystem=c("LPS","RAS"),labels=dataname) {
    dataname <- deparse(substitute(x))
    if (!grepl("*.json$", filename)) 
        filename <- paste0(filename,".mrk.json")
    nrx <- nrow(x)
    locked <- TRUE
    type <- match.arg(type[1],c("Curve","Fiducial","ClosedCurve"))

    if (labels[[1]] == dataname)
        mylabels <- paste0(dataname,"-",1:nrx)
    else if (length(labels) == 1 || type != "Fiducial")
        mylabels <- paste0(labels,"-",1:nrx)
    else if (length(labels) == nrx)
        mylabels <- labels
    
       
    coordinateSystem <- match.arg(coordinateSystem[1],c("LPS","RAS"))
    
    ## setup markups
    orientation <- c(-1,0,0,0,-1,0,0,0,1)
    if (coordinateSystem == "RAS")
        orientation <- c(1,0,0,0,1,0,0,0,1)
    position <- lapply(1:nrx, function(y) y <- x[y,] )
    cp <- data.frame(id=as.character(1:nrx),label=mylabels,description=rep("",nrx),associatedNodeID=rep("",nrx))
    cp$position <- position
    cp$orientation <- lapply(1:nrx,function(x) x <- orientation)
    cp$selected <- rep(TRUE,nrx)
    cp$locked <- FALSE
    cp$visibility <- rep(TRUE,nrx)
    cp$positionStatus <- rep("defined",nrx)

    markups <- data.frame(type=type,coordinateSystem=coordinateSystem,locked=TRUE,labelFormat="%N-%d")
    markups$controlPoints=list(cp)
    markups$display=createDisplayInfo()

    
    out <- list("@schema"="https://raw.githubusercontent.com/slicer/slicer/master/Modules/Loadable/Markups/Resources/Schema/markups-schema-v1.0.0.json#")
    
    out$markups <- markups
    write_json(out,pretty=T,auto_unbox = T,filename,digits=NA,always_decimal=TRUE)
    
}


createDisplayInfo <- function(visibility=TRUE,opacity=1,color=c(0,4,1,1),selectedColor=c(1.0000000, 0.5000076 ,0.5000076),
                              propertiesLabelVisibility=TRUE,pointLabelsVisibility=FALSE,textScale=3,glyphType="Sphere3D",
                              sliceProjectionUseFiducialColor=TRUE,sliceProjectionOutlinedBehindSlicePlane=FALSE,
                              glyphScale=1, glyphSize=5, useGlyphScale=TRUE, sliceProjection = FALSE,
                              sliceProjectionColor=c(1,1,1),sliceProjectionOpacity=0.6,lineThickness=0.2,
                              lineColorFadingStart=1,lineColorFadingEnd=10,lineColorFadingSaturation=1,
                              lineColorFadingHueOffset=0,handlesInteractive=FALSE,snapMode="toVisibleSurface") {


    display <- data.frame(visibility=visibility,opacity=opacity)
    display$color <- list(color)
    display$selectedColor <- list(selectedColor)
    display2 <- data.frame(propertiesLabelVisibility=propertiesLabelVisibility,pointLabelsVisibility=pointLabelsVisibility,textScale=textScale,glyphType=glyphType,sliceProjectionUseFiducialColor=sliceProjectionUseFiducialColor,sliceProjectionOutlinedBehindSlicePlane=sliceProjectionOutlinedBehindSlicePlane,glyphScale=glyphScale, glyphSize=glyphSize, useGlyphScale=useGlyphScale, sliceProjection = sliceProjection)
    ## combine
    display <- cbind(display,display2)
    
    display$sliceProjectionColor <- list(c(1,1,1))

    display2 <- data.frame(sliceProjectionOpacity=sliceProjectionOpacity,lineThickness=lineThickness,lineColorFadingStart=lineColorFadingStart,lineColorFadingEnd=lineColorFadingEnd,lineColorFadingSaturation=lineColorFadingSaturation,lineColorFadingHueOffset=lineColorFadingHueOffset,handlesInteractive=handlesInteractive,snapMode=snapMode)
    display <- cbind(display,display2)

    return(display)
}
