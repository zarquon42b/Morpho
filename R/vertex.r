unrefVertex <- function(mesh)
  {
    it <- mesh$it
    vind <- 1:dim(mesh$vb)[2]
    unref <- which(! vind %in% it)
    return(unref)
  }

rmVertex <- function(mesh,index,keep=FALSE)
  {
      if (! keep)
          {
      it <- mesh$it
     
      itdim <- dim(it)
      lRm <- length(index)
   
      vbn <- dim(mesh$vb)[2]
      indOrig <- 1:vbn
      indOut <- indOrig*0
      indNew <- 1:(vbn-lRm)     
      indOut[-index] <- indNew
#print(indOut)
    
    facefun <- function(x)
      {
        x <- indOut[x]
        return(x)
      }
    if (!is.null(it))
      {
        it <- matrix(facefun(it),itdim)
        checkface <- rep(0,itdim[2])
        storage.mode(it) <- "integer"
        storage.mode(checkface) <- "integer"
        
        checkface <- .Fortran("face_zero",it,itdim[2],checkface)[[3]]
        invalface <- which(checkface == 0) #;print(invalface)
        if (length(invalface) > 0)
          {
            mesh$it <- it[,-invalface]
            if(!is.null(mesh$material$color))
              mesh$material$color <- mesh$material$color[,-invalface]
          }
        else
          {
            mesh$it <- it
          }
        if (0 %in% dim(it))
          {mesh$it <- NULL
         }
      }
      
      mesh$vb <- mesh$vb[,-index]
      mesh <- adnormals(mesh)
  }
      else
          mesh <- rmVertex(mesh,c(1:ncol(mesh$vb))[-sort(index)])
      return(mesh)

  }
vert2points <- function(mesh)
  {
    out <- t(mesh$vb[1:3,])
    return(out)
  }
rmUnrefVertex <- function(mesh)
  {
    unref <- unrefVertex(mesh)
    lunr <- length(unref)
    
    if (lunr > 0)
      {
        mesh <- rmVertex(mesh,unref)
        cat(paste(" removed",lunr,"unreferenced vertices\n"))
      }
    else
      {
        cat(" nothing to remove\n")
      }
    return(mesh)
  }
