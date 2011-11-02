file2mesh<-function(file,clean=TRUE)
 
{
  noclean=NULL
  if (!clean)
    {noclean <- " --noclean"
   }
  system(paste("ply2ascii ",file," dump.ascii.ply",noclean,sep=""))
  mesh<-ply2mesh("dump.ascii.ply")
  unlink("dump.ascii.ply")
  return(mesh)
}
	
