find.faces <- function(mesh,a,b,c)
  {
    it <- mesh$it
    b1 <- apply(it,2,function(x){x <- a%in% x ;return(x)})
    b2 <- apply(it,2,function(x){x <- b%in% x ;return(x)})
    b3 <- apply(it,2,function(x){x <- c%in% x ;return(x)})
    out <- which((b1*b2*b3)==1)
    return(out)
  }
del.faces <- function(mesh,faces)
  {
    flen <- length(faces)
    mat <- matrix(faces,3,flen/3)
    count <- 1
    del <- NULL
    for (i in 1:dim(mat)[2])
      {tmp <- find.faces(mesh,mat[1,i],mat[2,i],mat[3,i])
       if (length(tmp == 1))
          {
            del[count] <- tmp
            count <- count+1
          }
     }
    print(del)
    mesh$it <- mesh$it[,-del]
    return(mesh)
  }
add.faces <- function(mesh,faces)
  {
    flen <- length(faces)
    mat <- matrix(faces,3,flen/3)
    mesh$it <- cbind(mesh$it,mat)
    return(mesh)
  }
