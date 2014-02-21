addo <- function(arr) {
    if (!is.numeric(arr) || length(dim(arr)) != 3)
        stop("please provide 3D numeric array")
    out <- .Call("addo",arr)
    return(out)
}
arrMean3 <- function(arr) {
    if (!is.numeric(arr) || length(dim(arr)) != 3)
        stop("please provide 3D numeric array")
    out <- .Call("arrMean3",arr)
    return(out)
}
