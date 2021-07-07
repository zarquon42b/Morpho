#' update a vector of indices after removal of some referenced items 
#'
#' update a vector of indices after removal of some referenced items
#' @param x vector containing indices (e.g. to matrix rows)
#' @param indexrange maximum range of the index in the referenced item structure
#' @param ignore integer vector: remove those items from the original structure
#'
#' @examples
#' refItem <- matrix(1:10,5,2)
#' index <- c(1,3,5) # this indexes some rows of the matrix we are interested in
#' ## now we want to ignore row 2 and 5 and want to update the index so it will still fit
#' indexNew <- updateIndices(index,c(2,5),indexrange=5)
#'
#' ## Here a more useful example:
#' data(boneData)
#' left <- c(4,6,8)
#'   ## determine corresponding Landmarks on the right side:
#'     # important: keep same order
#'     right <- c(3,5,7)
#'     pairedLM <- cbind(left,right)
#' ## now we want to remove some landmarks and need to updated the pairedLM indices
#' ignore <- c(5,6)
#' mynewboneLM <- boneLM[-ignore,,]
#' pairedLMnew <- apply(pairedLM,2,updateIndices,ignore=ignore,indexrange=dim(boneLM)[1])
#' @export
updateIndices <- function(x,ignore,indexrange) {
    k <- indexrange
    li <- length(ignore)
    lm.old <- c(1:k)[-ignore]
    mat.ptr <- matrix(c(1:(k - li), lm.old), k - li, 2)
    ptr <- function(xo) {
        if (length(which(ignore %in% xo)) != 0) 
            xo <- xo[-which(xo %in% ignore)]
        for (i in 1:(k - li)) xo[which(xo == mat.ptr[i, 2])] <- mat.ptr[i,1]
        return(xo)
    }
    return(ptr(x))
    
} 
