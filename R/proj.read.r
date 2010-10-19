proj.read<-function(lm,mesh)
{	if (is.character(mesh))
		{proj.back(lm,mesh)
		}
	
	else 
		{mesh2ply(mesh,"dump0")
		proj.back(lm,"dump0.ply")
		unlink("dump0.ply")
		}
	
	data<-ply2mesh("out_cloud.ply")
	unlink("out_cloud.ply")
	return(data)
}
