warpmovie3d<-function(x,y,n,col=skin1,palindrome=FALSE,folder=NULL,movie="warpmovie",add=F,...)
{	wdold<-getwd()
	if(!is.null(folder))
		{
			if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/")
			{folder<-paste(folder,"/",sep="")
                         dir.create(folder,showWarnings=F)
			setwd(folder)
			}
		}
        if (!add)
          {
            open3d()
          }
                                        #	points3d(t(x$vb[1:3,]),col="white",cex=0)
	#print(folder)
	for (i in 0:n)
		{mesh<-x
		mesh$vb[1:3,]<-(i/n)*y$vb[1:3,]+(1-(i/n))*x$vb[1:3,]
		a <- shade3d(mesh,col=col,...)
		if (i ==0)
			{readline("please select view and press return\n")
			}
	
		filename <- sprintf("%s%03d.png", movie, i)
		rgl.snapshot(filename,fmt="png")
		rgl.pop("shapes",id=a)
		}
	
	if (palindrome) ## go the other way ##
		{
		for (i in 1:(n-1))
			{mesh<-x
			mesh$vb[1:3,]<-(i/n)*x$vb[1:3,]+(1-(i/n))*y$vb[1:3,]
			a <- shade3d(mesh,col=col,...)
			filename <- sprintf("%s%03d.png", movie, i+n)
			rgl.snapshot(filename,fmt="png")
			rgl.pop("shapes",id=a)	
			}
		}
        rgl.close()
	setwd(wdold)
}
