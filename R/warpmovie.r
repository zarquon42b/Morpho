warpmovie<-function(x,y,n,col=skin1,palindrome=FALSE,folder="",movie="warpmovie")
{	if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/")
		{folder<-paste(folder,"/",sep="")
		}
	wdold<-getwd()
	setwd(folder)
	#print(folder)
	for (i in 0:n)
		{mesh<-x
		mesh$vb[1:3,]<-(i/n)*y$vb[1:3,]+(1-(i/n))*x$vb[1:3,]
		shade3d(mesh,col=col)
		if (i ==0)
			{readline("please select view\n")
			}
	
		filename <- sprintf("%s%03d.png", movie, i)
		rgl.snapshot(filename,fmt="png")
		rgl.clear()
		}
	
	if (palindrome) ## go the other way ##
		{
		for (i in 1:(n-1))
			{mesh<-x
			mesh$vb[1:3,]<-(i/n)*x$vb[1:3,]+(1-(i/n))*y$vb[1:3,]
			shade3d(mesh,col=col)
			filename <- sprintf("%s%03d.png", movie, i+n)
			rgl.snapshot(filename,fmt="png")
			rgl.clear()	
			}
		}
	setwd(wdold)
}
