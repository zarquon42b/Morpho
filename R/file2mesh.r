#' @rdname ply2mesh
#' @importFrom Rvcg vcgImport
#' @export
file2mesh <- function(filename,clean=TRUE,readcol=FALSE) {
    mesh <- vcgImport(filename, clean=clean, readcol=readcol)
    return(mesh)
}
	
