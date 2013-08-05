warpmovie3d <- function (x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",...) UseMethod("warpmovie3d")
warpmovie3d.matrix <- function(x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",add=FALSE,close=TRUE,countbegin=0,ask=TRUE,radius=NULL,links=NULL,lwd=1,...)
{	#wdold <- getwd()
	if(!is.null(folder))
		{
			if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/")
			{folder <- paste(folder,"/",sep="")
                         dir.create(folder,showWarnings=F)
                         movie <- paste(folder,movie,sep="")
			#setwd(folder)
			}
		}
        if (!add)
          {
            open3d()
          }
         k <- dim(x)[1]
        if (is.null(radius))
          radius <- (cSize(x)/sqrt(k))*(1/80)
        ## get bbox
        bbox <- apply(rbind(x,y),2,range)
        bbox <- expand.grid(bbox[,1],bbox[,2],bbox[,3])
        points3d(bbox,col="white")
                            
	for (i in 0:n)
		{
                  mesh <- x
                  mesh  <- (i/n)*y+(1-(i/n))*x
                 
                  a <- spheres3d(mesh,col=col,radius=radius,...)
                  if (!is.null(links))
                    {a1 <- lineplot(mesh,links,col=col,lwd=lwd)
                    a <- append(a,a1)
                   }
                  if (i ==0 && ask==TRUE)
                    {readline("please select view and press return\n")
                   }
                  
                  filename <- sprintf("%s%04d.png", movie, countbegin+i)
                  rgl.snapshot(filename,fmt="png")
                  rgl.pop("shapes",id=a)
		}
	
	if (palindrome) ## go the other way ##
		{
		for (i in 1:(n-1))
			{mesh <- x
                         mesh <- (i/n)*x+(1-(i/n))*y
                         a <- spheres3d(mesh,col=col,radius=radius,...)
                         if (!is.null(links))
                           {a1 <- lineplot(mesh,links,col=col,lwd=lwd)
                            a <- append(a,a1)
                          }
                         filename <- sprintf("%s%04d.png", movie, countbegin+i+n)
                         rgl.snapshot(filename,fmt="png")
                         rgl.pop("shapes",id=a)	
                       }
		}
        if (close)
        rgl.close()
	#setwd(wdold)
}
warpmovie2d <- function(x,y,n,col="green",palindrome=FALSE,folder=NULL,movie="warpmovie",links=NULL,lwd=1,imagedim = "800x800",par=list(xaxt="n",yaxt="n",bty="n"),...)
{	wdold <- getwd()
        widxheight <- as.integer(strsplit(imagedim, split = "x")[[1]])
	if(!is.null(folder))
		{
			if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/")
			{folder <- paste(folder,"/",sep="")
                         dir.create(folder,showWarnings=F)
			setwd(folder)
			}
		}
       
         k <- dim(x)[1]
       
        ## get bbox
        bbox <- apply(rbind(x,y),2,range)
        bbox <- expand.grid(bbox[,1],bbox[,2])
       
       # plot(bbox,col="white")
                            
	for (i in 0:n)
		{
                  mesh <- x
                  mesh  <- (i/n)*y+(1-(i/n))*x
                  filename <- sprintf("%s%04d.png", movie, i)
                  png(filename,width = widxheight[1], height = widxheight[2])
                  par(par)
                  plot(bbox,asp=1,cex=0,xlab="",ylab="")
                  points(mesh,col=col,...)
                  if (!is.null(links))
                    lineplot(mesh,links,col=col,lwd=lwd)

                  dev.off()
                  
		}
	
	if (palindrome) ## go the other way ##
		{
		for (i in 1:(n-1))
                  {mesh <- x
                   mesh <- (i/n)*x+(1-(i/n))*y
                   filename <- sprintf("%s%03d.png", movie, i+n)
                   png(filename,width = widxheight[1], height = widxheight[2])
                   par(par)
                   plot(bbox,asp=1,cex=0)
                   points(mesh,col=col,...)
                  
                   if (!is.null(links))
                     lineplot(mesh,links,col=col,lwd=lwd)

                   dev.off()
                 }
              }
       	setwd(wdold)
}
