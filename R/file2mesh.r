file2mesh<-function(file)
{      system(paste("ply2ascii ",file," dump.ascii.ply",sep=""))
	mesh<-ply2mesh("dump.ascii.ply")
	unlink("dump.ascii.ply")
	return(mesh)
}
	
