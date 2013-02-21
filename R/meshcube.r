meshcube <- function(x)
  {
    bbox <- apply(vert2points(x), 2, range)
    bbox <- expand.grid(bbox[, 1], bbox[, 2], bbox[, 3])
    return(bbox)
  }
