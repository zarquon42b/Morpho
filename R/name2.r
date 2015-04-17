#' extract data from array names
#' 
#' extract data from array names and convert to factors or numbers
#' 
#' If an array is used as input, the data info is expected to be in the 3rd
#' dimension, for a matrix, rownames are used.
#' 
#' @title extract data from array names
#' @param x data, can be a three-dimensional array, a matrix, a named list or a
#' vector containing names to split
#' @param sep character by which to split the strings
#' @param which integer or vector of integers, if more entries are selected,
#' they will be concatenated by the string specified with the option
#' 'collapse'.
#' @param collapse character by which to collapse data if two strings are to be
#' concatenated
#' @param dif logical: calculate difference if two fields containing numbers
#' are selected.
#' @param as.factor logical: if TRUE, a factor vector will be returned, strings otherwise.
#' @return returns a vector containing factors or numbers
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' 
#' data <- matrix(rnorm(200),100,2)
#' id <- paste("id",1:100,sep="")
#' pop <- c(rep("pop1",50),rep("pop2",50))
#' sex <- c(rep("male",50),rep("female",50))
#' age <- floor(rnorm(100,mean=50,sd=10))
#' rownames(data) <- paste(id,pop,sex,age,sep="_")
#' infos <- data.frame(pop=name2factor(data,which=2))
#' infos$age <- name2num(data,which=4)
#' infos$pop.sex <- name2factor(data,which=2:3)
#' 
#' 
#' @export
name2factor <- function(x,sep="_",which,collapse=sep, as.factor=TRUE) {
    if (length(dim(x))==3) {
        names <- dimnames(x)[[3]]
    } else if (length(dim(x))==2) {
          names <- dimnames(x)[[1]]
      } else if (is.list(x)) {        
            names <- names(x)
        } else {
              names <- x
          }
    fac <- strsplit(names,split=sep)
    for (i in 1:length(fac)) {                 
        fac[[i]] <- paste(fac[[i]][which],collapse=collapse)
    }
    
    fac <- unlist(fac)
    if (as.factor)
        fac <- factor(fac)
    return(fac)
}
#' @rdname name2factor
#' @export
name2num <- function(x,sep="_",which,collapse=sep,dif=TRUE) {
    if (length(dim(x))==3) {
        names <- dimnames(x)[[3]]
    } else if (length(dim(x))==2) {
          names <- dimnames(x)[[1]]
      } else if (is.list(x)) {        
            names <- names(x)
        } else {
              names <- x
          }
    fac <- strsplit(names,split=sep)
    
    for (i in 1:length(fac)) {
        if (dif && length(which) > 1) {
            fac[[i]] <- diff(as.numeric(fac[[i]][which]))
        } else {
              fac[[i]] <- as.numeric(fac[[i]][which])
          }
    }
    fac <- (unlist(fac))
    return(fac)
}
