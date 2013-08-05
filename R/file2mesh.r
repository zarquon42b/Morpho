file2mesh <- function(filename,clean=TRUE,readcol=FALSE)
 
{
  options=NULL
  if (!clean)
    {options <- paste(options,"--noclean")
   }
  if (readcol)
    {options <- paste(options,"--color")
   }
  system(paste("ply2ascii ",filename," dump.ascii.ply",options,sep=""))
  mesh <- ply2mesh("dump.ascii.ply",readcol=readcol)
  unlink("dump.ascii.ply")
  return(mesh)
}
	
