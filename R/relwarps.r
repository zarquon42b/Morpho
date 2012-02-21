relWarps<-function(data,scale=TRUE,CSinit=TRUE,alpha=1,tol=1e-10,orp=FALSE)
{	n<-dim(data)[3]
	m<-dim(data)[2]
	k<-dim(data)[1]

	### superimpose data ###
	proc<-sc.procGPA(data,scale=scale,CSinit=CSinit)
	if (orp)
          {
            proc$rotated<-orp(proc$rotated)
          }
        proc$mshape <- apply(proc$rotated,1:2,mean)
	### create bending energy matrix ###
	if (m==2)
		{BE<-CreateL2D(proc$mshape)$Lsubk
		}
	else
		{BE<-(CreateL(proc$mshape)$Lsubk)
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
	eigBE<-eigen(BE,symmetric=TRUE)
	
	eigBE$values<-Re(eigBE$values)
	eigBE$vectors<-Re(eigBE$vectors)
	
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
		eigCOVCOM<-eigen(covcom,symmetric=TRUE)
		nonz<-which(eigCOVCOM$values>tol)
		bescores<-t(t(eigCOVCOM$vectors[,nonz])%*%BE2%*%t(vecs))

	
	### calculate uniform component scores ###
		U<-NULL
		uniscores<-NULL
        
	if (m==1)
		{
                  rotms<-eigen(crossprod(proc$mshape))$vectors
                  if (det(rotms) < 0)
                    {rotms[,1]<-rotms[,1]*-1
                                        #rotms<-rotms*c(1,-1)#%*%matrix(c(0,-1,1,0),2,2)
                   }
		print(dim(rotms))
		msrot<-(proc$mshape%*%rotms)/cSize(proc$mshape)
		al<-sum(msrot[,1]^2)
		be<-sum(msrot[,2]^2)
		ga<-sqrt(al*be)
		u1<-c(al*msrot[,2],be*msrot[,1])/ga
		u2<-c(-be*msrot[,1],al*msrot[,2])/ga
                  U<-cbind(u1,u2)
                 
                  uniscores<-vecs%*%U
		}
        else
          {
            vecData <- vecx(proc$rotated)
### Rohlf first method ###
            
            E <- eigBE$vectors[,-zero]
            N <- diag(rep(1,k))-E%*%solve(crossprod(E))%*%t(E)
            V <- t(t(vecData)-as.vector(proc$mshape))
          #  NIk <- matrix(0,k+m,k+m)
            NIk <- kronecker(E,diag(rep(1,m)))
            
         #   for (i in 1:m)
         #     {
               
         #       NIk[((i-1)*k+1):(i*k),((i-1)*k+1):(i*k)] <- N
         #     }
            print(c(dim(V),dim(N)))
            svdBend <- svd(V%*%NIk)
            LS <- svdBend$u%*%diag(svdBend$d)
            uniscores <- LS[,1:(m+0.5*m*(m-1)-1)]

####Rohlf second method
          
           # Bxyz <- NULL
           # tXcInv <- solve(crossprod(proc$mshape))
           # for(i in 1:m)
           #   {  
           #     Bxyz <- cbind(Bxyz,t(tXcInv%*%t(proc$mshape)%*%t(vecData[,((i-1)*k+1):(i*k)])))
           #     print(dim(Bxyz))
           #   }
           # Bend <- cbind(t(Bxyz[[1]]),t(Bxyz[[2]]),t(Bxyz[[3]]))
           # print(dim(Bxyz))
           # svdBend <- svd(Bxyz)
           # LS <- svdBend$u%*%diag(svdBend$d)
           # print(dim(LS))
           # uniscores <- LS[,1:(m+0.5*m*(m-1)-1)]
          }
        

### create Variance table according to eigenvalues ###
		values<-eigCOVCOM$values[nonz]
		if (length(values)==1)
          	{Var<-values}
        else
        	{
          	Var<-matrix(NA,length(values),3)
          	Var[,1]<-values
        
          	for (i in 1:length(values))
            		{
              		Var[i,2]<-(values[i]/sum(values))*100
            		}
          	Var[1,3]<- Var[1,2]
          	for (i in 2:length(values))
           		{         
             		Var[i,3]<-Var[i,2]+ Var[i-1,3]
            		}
          	colnames(Var)<-c("eigenvalues","% Variance","Cumulative %")
        	}
      
	
	return(list(bescores=bescores,uniscores=uniscores,BE2=BE2,vecs=vecs,U=U,covSc=eigCOVCOM,invBE2=invBE2,proc=proc,Var=Var))
}
	
	

	
