scaleproc <- function (arr) {
    if (!is.numeric(arr) || length(dim(arr)) != 3)
        stop("please provide 3D numeric array")
    out <- .Call("scaleprocCpp",arr)
    return(out)
}

