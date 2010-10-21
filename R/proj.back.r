proj.back<-function(data,surface,dataname=NULL,outname=NULL)
{	
	if (is.null(dataname))
		{dataname<-"out"}
	write.obj(cbind("v",data),filename=dataname)
	if (is.null(outname))
		{command<-paste("trimesh_project"," ",dataname,".obj"," ",surface,sep="")
		}
	else
		{command<-paste("trimesh_project"," ",dataname,".obj"," ",surface," ",outname,sep="")
		}
	
	#command<-paste("trimesh_project"," ",dataname," ",surface," ",outname,sep="")	
	system(command)
	unlink(paste(dataname,".obj",sep="")) #clean up
	
}
