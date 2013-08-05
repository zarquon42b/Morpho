scalemesh <- function(mesh,size,center=c("bbox","mean", "none"))
{	
    getmean <- TRUE
  if (substr(center[1],1L,1L) =="b")
    meshmean <- apply(meshcube(mesh),2,mean)
  else if (substr(center[1],1L,1L) =="m")
    meshmean <- apply(vert2points(mesh),2,mean)
  else if (substr(center[1],1L,1L) =="n")
     getmean <- FALSE
  else
    stop("Please provide valid centering method\n")
    if (getmean)
        {
            mesh <- translate3d(mesh,-meshmean[1],-meshmean[2],-meshmean[3])
            mesh$vb[1:3,] <- mesh$vb[1:3,]*size
            mesh <- translate3d(mesh,meshmean[1],meshmean[2],meshmean[3])
        }
   else
       mesh$vb[1:3,] <- mesh$vb[1:3,]*size
  return(mesh)
}
	
