proj.back<-function(data,surface)
{	write.obj(cbind("v",data),filename="out")
	command<-paste("trimesh_project"," ","out.obj"," ",surface,sep="")
	system(command)
	unlink("out.obj") #clean up
	
}
