warp_obj<-function(obj,matr,matt,save=TRUE,filename="default")
{
      
      vertfr<-subset(obj,obj[,1]=="v")
      trifr<-subset(obj,obj[,1]=="f")
      vert<-as.matrix(subset(obj,obj[,1]=="v")[,2:4])
      tri<-as.matrix(subset(obj,obj[,1]=="f")[,2:4])
      
      warp<-tps3d(vert,matr,matt)
      obj[which(obj[,1]=="v"),][,2:4]<-warp
      
      
      
      sv1<-svd(t(matt)%*%matr)
      if (sign(det(sv1$v))<0)
      {obj<-conv2backf(obj)
        cat("reflection is involved")}
      
      if (save==TRUE)
        {write.obj(obj,filename=filename)}
      
      return(obj)
}
      