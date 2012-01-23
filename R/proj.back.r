proj.back<-function(data,surface,dataname=NULL,outname=NULL,smooth=TRUE,ignore.stdout=FALSE)
{	
	smoothopt<-NULL
	if (!smooth)
	{smoothopt<-" --nosmooth "
	}
	if (is.null(dataname))
		{dataname<-"out"}
	write.obj(cbind("v",data),filename=dataname)
	if (is.null(outname))
		{command<-paste("trimesh_project"," ",dataname,".obj"," ",surface,smoothopt,sep="")
		}
	else
		{command<-paste("trimesh_project"," ",dataname,".obj"," ",surface," -o ",outname,smoothopt,sep="")
		}
	
	#command<-paste("trimesh_project"," ",dataname," ",surface," ",outname,sep="")	
	system(command,ignore.stdout=ignore.stdout)
	unlink(paste(dataname,".obj",sep="")) #clean up
	
}
