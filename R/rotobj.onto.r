rotobj.onto<-function(obj,refmat,tarmat)
{ 	rot<-rotonto(tarmat,refmat)
	centmesh<-t(apply(obj[which(obj[,1]=="v"),2:4],1,function(x){x - rot$transy}))
	centmeshr<-centmesh%*%rot$gamm
  	
	obj[which(obj[,1]=="v"),2:4]<-t(apply(centmeshr,1,function(x){x+rot$trans}))
  	
	if (sign(det(rot$gamm)<0))
  		{obj<-conv2backf(obj)
		}

  	return(list(obj=obj,yrot=rot$yrot))
}
