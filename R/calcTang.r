calcTang<-function(datamatrix,SMvector,outlines,deselect=FALSE)
{     if (is.list(outlines)==FALSE)
      {outlines<-list(outlines)}
      
      dims<-dim(datamatrix)[2]
      k<-dim(datamatrix)[1]
      m<-length(SMvector)
      tanvec<-matrix(0,k,dims)
      
      if (deselect==TRUE)
      {SMvector<-c(1:k)[-SMvector]}
      
      for ( j in 1:length(outlines))
        lt<-length(outlines[[j]])
        temp<-outlines[[j]]
        for (i in 1:lt)
          
          {if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt)
           {tanvec[temp[i],]<-(datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
           }
          
           
          else if (temp[i]%in%SMvector==TRUE && i==1)
           {tanvec[temp[i],]<-(datamatrix[temp[i],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i],]-datamatrix[temp[i+1],])^2))
           }
          else if (temp[i]%in%SMvector==TRUE && i==lt)
           {tanvec[temp[i],]<-(datamatrix[temp[i-1],]-datamatrix[temp[i],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i],])^2))
          }
          #else {tanvec[i,]<-c(rep(0,dims))}
          }
    return(list(tanvec,SMvector))
}