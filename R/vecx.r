vecx <- function(x,byrow=FALSE)
{  dims <- dim(x)
	n <- dims[3]
	k <- dims[1]
	m <- dims[2]
   names <- dimnames(x)[[3]]
   
	vecs <- matrix(0,n,k*m)
	for(i in 1:n)
          {
            if (byrow)
              {
                vecs[i,] <- as.vector(t(x[,,i]))
              }
            else
              {vecs[i,] <- as.vector(x[,,i])
             }
          }
   if (!is.null(names))
     {
       rownames(vecs) <- names
     }
	#vecs <- apply(vecs,2,scale,scale=F)
	return(vecs)
}
	
