#' calculate the corners of a mesh's bouning box
#' 
#' calculate the corners of a mesh's bouning box
#' 
#' 
#' @param x object of class 'mesh3d'
#' @return returns a 8 x 3 matrix with the coordinates of the corners of the
#' bounding box.
#' 
#' @examples
#' 
#' require(rgl)
#' data(boneData)
#' mc <- meshcube(skull_0144_ch_fe.mesh)
#' \dontrun{
#' spheres3d(mc)
#' wire3d(skull_0144_ch_fe.mesh)
#' }
#' 
#' @export
meshcube <- function(x)
  {
    bbox <- apply(vert2points(x), 2, range)
    bbox <- expand.grid(bbox[, 1], bbox[, 2], bbox[, 3])
    return(bbox)
  }
