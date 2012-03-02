relWarps<-function(data,scale=TRUE,CSinit=TRUE,alpha=1,tol=1e-10,orp=TRUE)
{
  n<-dim(data)[3]
  m<-dim(data)[2]
  k<-dim(data)[1]
  datanames <- dimnames(data)[[3]]
### superimpose data ###
  proc<-ProcGPA(data,scale=scale,CSinit=CSinit,silent=TRUE)
  if (orp)
    {
      proc$rotated<-orp(proc$rotated)
      proc$mshape <- apply(proc$rotated,1:2,mean)
    }
  
  
### create bending energy matrix ###
  if (m==2)
    {
      BE<-CreateL2D(proc$mshape)$Lsubk
    }
  else
    {
      BE<-(CreateL(proc$mshape)$Lsubk)
    }
### vectorize and scale superimposed data ###
  vecs <- vecx(proc$rotated)
  vecs<-apply(vecs,2,scale,scale=F)
  dimnames(proc$rotated)[[3]] <- rownames(vecs) <- datanames
  
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
    {
      BE2<-diag(rep(1,k*m))
      invBE2<-BE2
    }
  
### generate covariance structure of scaled space ###
  covcom<-BE2%*%Sc%*%BE2	
  eigCOVCOM<-eigen(covcom,symmetric=TRUE)
  nonz<-which(eigCOVCOM$values>tol)
  bescores<-t(t(eigCOVCOM$vectors[,nonz])%*%BE2%*%t(vecs))
  rownames(bescores) <- rownames(vecs)
  bePCs <- IM %x% eigBE$vectors
  bePCs <-bePCs %*% diag(rep(diaginv,3)) %*% t(bePCs) %*%  eigCOVCOM$vectors
  
### calculate uniform component scores ###
  U<-NULL
  uniscores<-NULL
  
### Rohlf first method ###
  
  E <- eigBE$vectors[,-zero]
  N <- diag(rep(1,k))-E%*%solve(crossprod(E))%*%t(E)
  V <-vecs
  NIk <- IM %x% N
   
  svdBend <- svd(V%*%NIk)
  LS <- svdBend$u%*%diag(svdBend$d)
  uniscores <- LS[,1:(m+0.5*m*(m-1)-1)]
  rownames(uniscores) <- datanames   
  
### create Variance table according to eigenvalues ###
  values<-eigCOVCOM$values[nonz]
  if (length(values)==1)
    {
      Var<-values
    }
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
  
  
  return(list(bescores=bescores,uniscores=uniscores,Var=Var,mshape=proc$mshape,rotated=proc$rotated,bePCs=bePCs,uniPCs=svdBend$v[,1:(m+0.5*m*(m-1)-1)]))
  
}

