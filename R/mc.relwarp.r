mc.relwarps<-function(data,scale=TRUE,CSinit=TRUE,alpha=1,tol=1e-10)
{	n<-dim(data)[3]
	m<-dim(data)[2]
	k<-dim(data)[1]

	### superimpose data ###
	proc<-mc.procGPA(data,scale=scale,CSinit=CSinit)
	
	### create bending energy matrix ###
	if (m==2)
		{BE<-CreateL2D(proc$mshape)$Lsubk
		}
	else
		{BE<-CreateL(proc$mshape)$Lsubk
		}
	### vectorize and scale superimposed data ###
	vecs<-matrix(0,n,k*m)
	for(i in 1:n)
		{vecs[i,]<-as.vector(proc$rotated[,,i])
		}
	vecs<-apply(vecs,2,scale,scale=F)
	
	### generate covariance matrix of superimposed data ###
	Sc<-cov(vecs)
	
	### explore eigenstructure of BE ###	
	eigBE<-eigen(BE)
	zero<-which(eigBE$values<tol)
	diaginv<-diagBE<-eigBE$values*0
	diagBE[-zero]<-eigBE$values[-zero]^(-alpha/2)
	diaginv[-zero]<-eigBE$values[-zero]^(alpha/2)
	IM<-diag(rep(1,m))
	
	if (alpha !=0)
		{	
		BE2<-IM%x%(eigBE$vectors%*%diag(diagBE)%*%t(eigBE$vectors))
	
		invBE2<-IM%x%(eigBE$vectors%*%diag(diaginv)%*%t(eigBE$vectors))
		}
	else
		{BE2<-diag(rep(1,k*m))
		invBE2<-BE2
		}
	
	### generate covariance structure of scaled space ###
		covcom<-BE2%*%Sc%*%BE2	
		eigCOVCOM<-eigen(covcom)
		nonz<-which(eigCOVCOM$values>tol)
	
		#finalmat<-t(eigCOVCOM$vectors[,nonz])%*%BE2
		bescores<-t(t(eigCOVCOM$vectors[,nonz])%*%BE2%*%t(vecs))

	
	### calculate uniform component scores ###
		U<-NULL
		uniscores<-NULL
	if (m==2)
		{rotms<-eigen(crossprod(proc$mshape))$vectors
		if (det(rotms) < 0)
			{rotms[,1]<-rotms[,1]*-1
			#rotms<-rotms*c(1,-1)#%*%matrix(c(0,-1,1,0),2,2)
			}
		
		msrot<-(proc$mshape%*%rotms)/c.size(proc$mshape)
print(t(msrot[,1])%*%msrot[,2])
		plot(msrot,asp=1)
		al<-sum(msrot[,1]^2)
		be<-sum(msrot[,2]^2)
		ga<-sqrt(al*be)
		u1<-c(al*msrot[,2],be*msrot[,1])/ga
		u2<-c(-be*msrot[,1],al*msrot[,2])/ga
		U<-cbind(u1,u2)
		uniscores<-vecs%*%U
		
		
		}
	
	return(list(bescores=bescores,uniscores=uniscores,BE2=BE2,vecs=vecs,U=U,covSc=eigCOVCOM,invBE2=invBE2,proc=proc,rotms=rotms))
}
	
	

	
