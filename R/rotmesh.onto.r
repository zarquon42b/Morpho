rotmesh.onto<-function(mesh,refmat,tarmat,adnormals=TRUE,scale=FALSE)
{ 	rot<-rotonto(tarmat,refmat,scale=scale)
	
	transymat<-diag(c(rep(1,4)))
	transymat[4,1:3]<--rot$transy
	transmat<-diag(c(rep(1,4)))
	transmat[4,1:3]<-rot$trans
	#print(transmat)
	mesh$vb<-(apply(t(mesh$vb),1,function(x){x%*%transymat}))
	#print(mesh$vb[,1:2])	
	mesh$vb[1:3,]<-t(t(mesh$vb[1:3,])%*%rot$gamm)
	if (scale)
		{mesh$vb[1:3,]<-mesh$vb[1:3,]*rot$bet
		}
  	#mesh$vb[1:3,]<-t(centmeshr)
	mesh$vb<-(apply(t(mesh$vb),1,function(x){x%*%transmat}))
  	
	if (sign(det(rot$gamm)<0))
  		{mesh<-conv2backf(mesh)
		}
	if (adnormals)
		{mesh<-adnormals(mesh)
		}

  	return(list(mesh=mesh,yrot=rot$yrot))
}
