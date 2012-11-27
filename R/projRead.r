projBack<-function(data,surface,dataname=NULL,outname=NULL,smooth=TRUE,ignore.stdout=FALSE,sign=FALSE)
{	
	options<-NULL
        
	if (!smooth)
	{options <- paste(options,"--nosmooth")
	}
        if (sign)
          {
            options <- paste(options,"--sign")
          }
	if (is.null(dataname))
		{dataname<-"out"}
	write.obj(cbind("v",data),filename=dataname)
	if (is.null(outname))
		{command<-paste("trimesh_project"," ",dataname,".obj"," ",surface,options,sep="")
		}
	else
		{command<-paste("trimesh_project"," ",dataname,".obj"," ",surface," -o ",outname,options,sep="")
		}
	#print(command)
	#command<-paste("trimesh_project"," ",dataname," ",surface," ",outname,sep="")	
	system(command,ignore.stdout=ignore.stdout)
	unlink(paste(dataname,".obj",sep="")) #clean up
	
}
projRead<-function(lm,mesh,readnormals=TRUE,clean=TRUE,smooth=TRUE,ignore.stdout=FALSE,sign=FALSE,outname=NULL)
{	if (is.character(mesh))
		{projBack(lm,mesh,ignore.stdout=ignore.stdout,sign=sign,smooth=smooth,outname=outname)
		}
	
	else 
		{mesh2ply(mesh,"dump0")
		projBack(lm,"dump0.ply",smooth=smooth,ignore.stdout=ignore.stdout,sign=sign,outname=outname)
		unlink("dump0.ply")
		}
	
	data<-ply2mesh("out_cloud.ply",readnormals=readnormals)
	if (clean)
	{	
	unlink("out_cloud.ply")
	}
	return(data)
}
