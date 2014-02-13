#' estimate missing landmarks from their bilateral counterparts
#'
#' estimate missing landmarks from their bilateral counterparts
#'
#' @param x a matrix or an array containing landmarks (3D or 2D)
#' @param pairedLM  a k x 2 matrix containing the indices (rownumbers) of the
#' paired LM. E.g. the left column contains the lefthand landmarks, while the
#' right side contains the corresponding right hand landmarks.
#' @note in case both landmarks of a bilateral pair are missing a message will be issued. As well if there are missing landmarks on the midsaggital plane are detected.
#' @details the configurations are mirrored and the relabled version is matched onto the original using a thin-plate spline deformation. The missing landmark is now estimated using its bilateral counterpart.
#' @return a matrix or array with fixed missing bilateral landmarks.
#' @examples
#' data(boneData)
#' left <- c(4,6,8)
#' ## determine corresponding Landmarks on the right side:
#' # important: keep same order
#' right <- c(3,5,7)
#' pairedLM <- cbind(left, right)
#' exampmat <- boneLM[,,1]
#' exampmat[4,] <- NA #set 4th landmark to be NA
#' fixed <- fixLMmirror(exampmat, pairedLM=pairedLM)
#' \dontrun{
#' deformGrid3d(fixed, boneLM[,,1],ngrid=0)
#' ## result is a bit off due to actual asymmetry
#' }
#' @export
fixLMmirror <- function(x, pairedLM) UseMethod("fixLMmirror")

#' @rdname fixLMmirror
#' @export
fixLMmirror.array <- function(x,pairedLM) {
    n <- dim(x)[3]
    out <- x
    for (i in 1:n) {
        out[,,i] <- fixLMmirror(x[,,i],pairedLM=pairedLM)
    }
    return(out)
}

#' @rdname fixLMmirror
#' @export
fixLMmirror.matrix <- function(x,pairedLM) {
    data <- x
    m <- dim(x)[2]
    count <- 0
    k <- nrow(x)
    checklist <- NA
    unilatNA <- NULL
    for (j in 1:k) {
        if (TRUE %in% is.na(data[j,])) {
            count <- count+1
            checklist[count] <- j
        }
    }
    affected <- affectCol <- goodPaired <- NULL
    for (i in 1:nrow(pairedLM)) {
        if (prod(is.na(x[pairedLM[i,],]))) {
            warning(paste("paired landmarks",pairedLM[i,1], "and" ,pairedLM[i,2] , ": one landmark of each side must be present"))
            unilatNA <- append(unilatNA,pairedLM[i,])
            
        }
        checkPaired <- is.na(x[pairedLM[i,],])
        if (TRUE %in% checkPaired) {
            affected <- append(affected,i)
            checkCol <- (apply(checkPaired,1,prod))
            affectCol <- append(affectCol,which(!as.logical(checkCol)))
            goodPaired <- append(goodPaired,pairedLM[i,which(!as.logical(checkCol))])
        }
    }
    if (is.null(affected)) {
        message("no missing bilateral landmarks")
        return(x)
    }
   
    if (!prod(checklist %in% pairedLM)){
        warning("missing landmarks are not bilateral")
        unilatNA <- append(unilatNA,checklist[which(! checklist %in% pairedLM)])
    }
        
    xmir <- x %*% diag(c(-1,1,1))[1:m, 1:m]#mirror landmarks
    xmir[c(pairedLM),] <- xmir[c(pairedLM[,2:1]),]##relabel landmarks
    xref <- xmir[-c(unilatNA,pairedLM[affected,]),]
    xtar <- x[-c(unilatNA,pairedLM[affected,]),]
    xout <- tps3d(xmir,xref,xtar)
    xout[goodPaired,] <- x[goodPaired,]
    return(xout)

}
