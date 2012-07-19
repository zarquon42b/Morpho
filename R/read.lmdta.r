read.lmdta<-function(file="x")
{         x<-file
          A<-readLines(x)
          em<-which(A=="")
          idnames<-A[c((em[1]+1):(em[2]-1))]
          info<-strsplit(A[3]," ")[[1]]
          n2<-nchar(info[2])-1
          nspeci<-as.numeric(substr(info[2],1L,n2))
          ndim<-as.numeric(substr(info[6],5,nchar(info[6])))
          nlms<-as.numeric(info[3])/ndim
          eot<-em[2]
          B<-as.matrix(read.table(x,skip=eot),na.strings=as.numeric(info[5]))
          tt<-array(t(B),dim=c(ndim,nlms,nspeci))
          arr<-array(NA,dim=c(nlms,ndim,nspeci))
          for (i in 1:nspeci){arr[,,i]<-t(tt[,,i])}
          dimnames(arr)[[3]]<-idnames
          return(list(arr=arr,info=info,idnames=idnames))
}
