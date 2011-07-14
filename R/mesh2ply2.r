mesh2ply2<-function(x,filename="default",curvature=TRUE,neighbour=2)
{ 
  if (is.matrix(x))
    {x<-list(vb=x)}
  filename<-paste(filename,".ply2",sep="")
  vert<-x$vb[1:3,]
  vert<-round(vert,digits=6)
  if (!is.null(x$it))
    {face<-x$it-1
     fd<-3
     fn<-dim(face)[2]
  }
  else 
    {fn<-0
   } 
  vert.all<-vert
  vn<-dim(vert)[2]
  vn.all<-3
  texfile<-x$TextureFile
 
 
  
  header <- c(vn,fn)
  if (curvature)
    {header <- c(vn,fn,neighbour,1)
     faceraw <- face
   }
  else
    {
      faceraw<-rbind(fd,face)
    }
### start writing to file ###
  
  
  
### write vertex infos to header ###
  
  
  write.table(header,file=filename,append=FALSE,sep="",quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
                                        # cat(paste(fn,"\n"),file=filename,append=TRUE,sep="")
  
### write face infos and texture infos to header and finish ### 
  
  
                                        # cat(fn,file=filename,append=TRUE) 
  
  
 
### write vertices ###
  write.table(format(as.vector(vert.all),
                     scientific=F,trim=T),file=filename,sep=" ",append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, na = "") 
  
  
### write face ### 
                                        # if (!is.null(x$it)){
  write.table(format(as.vector(faceraw),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")   
#  }
# else
                                        #  {write.table(format(t(rbind(fd,face)),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "") 
                                        #  }
  
}
