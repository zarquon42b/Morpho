plotCurves <- function(curve,text=FALSE,color="rainbow")
  {
    limits <- as.integer(levels((as.factor(curve$id))))
    if (color !="rainbow")
      {points3d(curve$vertices,col=color)
     }
    if (text || color == "rainbow")
      {
        for (i in limits)
          {
            if (color == "rainbow")
              {
                col <- i+1
                points3d(curve$vertices[which(curve$id==i),],col=col)
              }
            else
              {col <- color
             }
            if (text)
              {
                text3d(apply(curve$vertices[which(curve$id==i),],2,mean),texts=paste(i),col=col)
              }
          }
      }
  }

readCurves <- function(x)
  {
    tmp <- as.numeric(readLines(x,n=1))
    tmp1 <-  readLines(x,n=tmp+3)
    tmp2 <- strsplit(tmp1,split=" ")
    chckentr <- unlist(lapply(tmp2,length))
    vertex <- which(chckentr==4)
    lv <- length(vertex)
    vertices <- matrix(as.numeric(unlist(tmp2[vertex])),lv,4,byrow=TRUE)

    n <- as.integer(vertices[,4])
    vertices <- vertices[,1:3]

    return(list(vertices=vertices,id=n))
  }

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
  
### write face infos and texture infos to header and finish ### 

### write vertices ###
  write.table(format(as.vector(vert.all),scientific=F,trim=T),file=filename,sep=" ",append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, na = "") 
  
  
### write face ### 
                                       
  write.table(format(as.vector(faceraw),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")   
  
}

file2ply2 <- function(filein,fileout,curvature=TRUE,neighbour=2)
  {
    tmpmesh <- file2mesh(filein)
    mesh2ply2(tmpmesh,filename=fileout,curvature=curvature,neighbour=neighbour)
  }
