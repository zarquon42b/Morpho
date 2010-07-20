conv2backf<-function(obj)
{ 	
	if (!("mesh3d" %in% class(obj)))
		{obj[which(obj[,1]=="f"),][,2:4]<- obj[which(obj[,1]=="f"),][,c(4,3,2)]
   		write.obj(obj)
		}
	else 
		{obj$it<-obj$it[c(3,2,1),]
		}
  return(obj)}
