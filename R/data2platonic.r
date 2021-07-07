#' creates 3D shapes from data to be saved as triangular meshes
#' 
#' creates 3D shapes from 3-dimensional data that can be saved as triangular meshes
#' @param datamatrix k x 3 data matrix
#' @param shape a 3D shape
#' @param col color value
#' @param scale logical: whether to scale the data to unit sd.
#' @param scalefactor scale the resulting shapes.
#' @return returns all shapes merged into a single mesh
#' @examples
#' 
#' mymesh <- data2platonic(iris[iris$Species=="setosa",1:3],scalefactor=0.1)
#' mymesh <- mergeMeshes(mymesh,data2platonic(iris[iris$Species=="versicolor",1:3],
#'                       shape=Rvcg::vcgIcosahedron(),scalefactor=0.1,col="green"))
#' mymesh <- mergeMeshes(mymesh,data2platonic(iris[iris$Species=="virginica",1:3],
#'                       shape=Rvcg::vcgTetrahedron(),scalefactor=0.1,col="blue"))
#' \dontrun{
#' rgl::shade3d(mymesh)
#' ## save to disk
#' Rvcg::vcgPlyWrite(mymesh,filename="3D_Data.ply")
#' }
#' @export
data2platonic <- function(datamatrix,shape=Rvcg::vcgSphere(),col="red",scale=FALSE,scalefactor=1) {
    myplatonic <- shape
    myplatonic <- Morpho::scalemesh(myplatonic,center="none",size=scalefactor)
    if(scale)
        matrix <- scale(matrix)
    col1mesh <- rgb(t(col2rgb(col)), maxColorValue = 255)
    matmesh <- lapply(1:nrow(datamatrix), function(x) x <- myplatonic)
    matmesh <- lapply(1:nrow(datamatrix), function(x) x <- translate3d(matmesh[[x]], x = datamatrix[x, 1], y = datamatrix[x, 2], z = datamatrix[x,3]))
    matmesh <- Morpho::mergeMeshes(matmesh)
    matmesh$material$color <- rep(col1mesh, ncol(matmesh$vb))
    matmesh$normals <- NULL
    return(matmesh)
}
