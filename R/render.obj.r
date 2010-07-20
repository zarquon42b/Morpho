render.obj<-function(obj,wire="s",col=1,size=size,plot=T,shapeout=T)     #options: w=wireframe, s=shade, p=points
{     if (dim(obj)[2]==4)
       {   
          vert<-as.matrix(subset(obj,obj[,1]=="v")[,2:4])
          tri<-as.matrix(subset(obj,obj[,1]=="f")[,2:4])
      }
      
      if (dim(obj)[2]==5)
      { meshlist<-list(0)
        vert<-as.matrix(subset(obj,obj[,1]=="v")[,2:4])
        tri<-as.matrix(subset(obj,obj[,1]=="f")[,2:5])
        
        
      
      for (i in 1:dim(tri)[1])                                                                         
          { 
            if  (NA %in% tri[i,4]==T)
            {meshlist[[i]]<-tmesh3d(c(vert[tri[i,1],],vert[tri[i,2],],vert[tri[i,3],]),indices=c(1:3),homogeneous=F)}
            
            else
            {meshlist[[i]]<-qmesh3d(c(vert[tri[i,1],],vert[tri[i,2],],vert[tri[i,3],],vert[tri[i,4],]),indices=c(1:4),homogeneous=F)}
          }
      }
      
      else if  (dim(obj)[2]==4)
        {
        meshlist<-list(0)


        for (i in 1:dim(tri)[1])
          { meshlist[[i]]<-tmesh3d(c(vert[tri[i,1],],vert[tri[i,2],],vert[tri[i,3],]),indices=c(1:3),homogeneous=F)}
        }  
        shape<-shapelist3d(meshlist,plot=F)

        if (wire=="w"&& plot==TRUE)
            {wire3d(shape,color=col)}

        if (wire=="s" && plot==TRUE)
            {shade3d(shape,color=col)}
        if  (wire=="p" && plot==TRUE)
            {dot3d(shape,color=col)}


        if (shapeout==T){return(shape=shape)}

}