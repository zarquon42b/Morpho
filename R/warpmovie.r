warpmovie3d.mesh3d<-function(x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",add=F,close=TRUE,countbegin=0,ask=TRUE,...)
{	#wdold<-getwd()
  if(!is.null(folder))
    {
      if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/")
        {folder<-paste(folder,"/",sep="")
         dir.create(folder,showWarnings=F)
         movie <- paste(folder,movie,sep="")
                                        #setwd(folder)
       }
    }
  if (!add)
    {
      open3d()
    }
  ## get bbox
  bbox <- apply(rbind(vert2points(x),vert2points(y)),2,range)
  bbox <- expand.grid(bbox[,1],bbox[,2],bbox[,3])
  points3d(bbox,col="white",alpha=0)
  
  for (i in 0:n)
    {
      mesh<-x
      mesh$vb[1:3,]<-(i/n)*y$vb[1:3,]+(1-(i/n))*x$vb[1:3,]
      mesh <- adnormals(mesh)
      a <- shade3d(mesh,col=col,...)
      if (i ==0 && ask==TRUE)
        {readline("please select view and press return\n")
       }
      
      filename <- sprintf("%s%03d.png", movie, countbegin+i)
      rgl.snapshot(filename,fmt="png")
      rgl.pop("shapes",id=a)
    }
  
  if (palindrome) ## go the other way ##
    {
      for (i in 1:(n-1))
        {mesh<-x
         mesh$vb[1:3,]<-(i/n)*x$vb[1:3,]+(1-(i/n))*y$vb[1:3,]
         mesh <- adnormals(mesh)
         a <- shade3d(mesh,col=col,...)
         filename <- sprintf("%s%03d.png", movie, countbegin+i+n)
         rgl.snapshot(filename,fmt="png")
         rgl.pop("shapes",id=a)	
       }
    }
  if (close)
    rgl.close()
                                        #setwd(wdold)
}
