warpmovie3d <- function (x,y,n,col="white",palindrome=FALSE,folder=NULL,movie="warpmovie",add=F,...) UseMethod("warpmovie3d")
warpmovie3d.matrix<-function(x,y,n,col="white",palindrome=FALSE,folder=NULL,movie="warpmovie",add=F,radius=NULL,...)
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
         k<-dim(x)[1]
        if (is.null(radius))
          radius <- (cSize(x)/sqrt(k))*(1/80)
        ## get bbox
        bbox <- apply(rbind(x,y),2,range)
        bbox <- expand.grid(bbox[,1],bbox[,2],bbox[,3])
        points3d(bbox,col="white")
                            
	for (i in 0:n)
		{
                  mesh<-x
                  mesh <-(i/n)*y+(1-(i/n))*x
                 
                  a <- spheres3d(mesh,col=col,radius=radius,...)
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
                         mesh <- (i/n)*x+(1-(i/n))*y
                         a <- spheres3d(mesh,col=col,radius=radius,...)
                         filename <- sprintf("%s%03d.png", movie, i+n)
                         rgl.snapshot(filename,fmt="png")
                         rgl.pop("shapes",id=a)	
                       }
		}
        rgl.close()
	setwd(wdold)
}
