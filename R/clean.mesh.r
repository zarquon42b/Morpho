clean.mesh<-function(dat.array=0,sur.path="sur",sur.name=NULL,sur.type="ply")
{	
	write(paste("<!DOCTYPE FilterScript>\n","<FilterScript>\n"," <filter name=\"Remove Duplicated Vertex\"/>\n","</FilterScript>",sep=""),file="clean.mlx")
	if(length(sur.name)==0)
	{	
		sur.name<-dimnames(dat.array)[[3]]
		sur.name<-paste(sur.path,"/",sur.name,".",sur.type,sep="")
      	}
	cat(sur.name)	
	n<-length(sur.name)
	for (i in 1:n)
		{
			command<-paste("meshlabserver -i"," ",sur.name[i]," -o ",sur.name[i]," -s clean.mlx",sep="")
			system(command)
		}
	unlink("clean.mlx")
	
}
