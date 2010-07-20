adnormals<- function(x) 
{	#### this is basically the function addNormals from the rgl-package with an option for checking degenerated faces ####
	

	veclen <- function(v) sqrt(sum(v^2))
	normalize <- function(v) v/veclen(v)
  v <- x$vb
  
  # Make sure v is homogeneous with unit w
  if (nrow(v) == 3) v <- rbind(v, 1)
  else v <- t( t(v)/v[4,] )
  
  normals <- v*0
  v <- v[1:3,]
  
  if (!is.null(x$it)) {
    it <- x$it
    for (i in 1:ncol(it)) {

	norm.raw<-(crossp( v[, it[1, i]] - v[, it[3, i]],v[, it[2, i]] - v[, it[1, i]]))
	if (veclen(norm.raw)==0) 
		{normal<-c(0,0,0)
		}
	    	
	else 
		{normal <- normalize(crossp( v[, it[1, i]] - v[, it[3, i]],v[, it[2, i]] - v[, it[1, i]]))
		}
	
      for (j in 1:3) {
		
	
      	if (sum(normals[1:3, it[j,i]]*normal) < 0)
      	  {normals[, it[j,i]] <- normals[, it[j,i]] + c(-normal, 1)}
		
		
      	else
          normals[, it[j,i]] <- normals[, it[j,i]] + c(normal, 1)
	
	#if (veclen(normals[1:3, it[j,i]])!=0) {normals[1:3, it[j,i]]<-normalize(normals[1:3, it[j,i]])}
      }
    }
  }
  
  if (!is.null(x$ib)) {
    it <- x$ib
    for (i in 1:ncol(it)) {
	
	norm.raw<-(crossp( v[, it[1, i]] - v[, it[3, i]],v[, it[2, i]] - v[, it[1, i]]))
	if (veclen(norm.raw)==0) 
		{normal<-c(0,0,0)
		}
	      	
	else 
		{normal <- normalize(crossp( v[, it[1, i]] - v[, it[3, i]],v[, it[2, i]] - v[, it[1, i]]))
		}
      
      for (j in 1:4) {
      	if (sum(normals[1:3, it[j,i]]*normal) < 0)
      	  normals[, it[j,i]] <- normals[, it[j,i]] + c(-normal, 1)
      	else
          normals[, it[j,i]] <- normals[, it[j,i]] + c(normal, 1)
      }

    }
  }
  normals <- t( t(normals)/normals[4,] )
  x$normals <- normals
  x
}


