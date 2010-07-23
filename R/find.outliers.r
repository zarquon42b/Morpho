find.outliers<-function(A,color=4,lwd=1,lcol=2)
{   	
	raw<-A
	n<-dim(A)[3]
	k<-dim(A)[1]
	m<-dim(A)[2]
	A<-procGPA(A,distances=FALSE,pcaoutput=FALSE)

###from here on the same as mc.find.outliers ###	
	rho<-0
	for(i in 1:n)
		{rho[i]<-angle.calc(A$rotated[,,i],A$mshape)$rho
		}
	A$rho<-rho
	
	
    	#rgl.clear()
    	#rgl.bg(color = "white")
   	disti<-data.frame(c(1:n),A$rho)
    	disti.sort<-disti[order(disti[,2],decreasing=T),]
    	colnames(disti.sort)[1]<-"# in array"
    	rownames(disti.sort)<-c(1:n)
    	outlier<-NULL
    	t1 <- 1
    	lo<-1
	if (m==3)
		{op<-open3d()	
		}  	
	
	while (t1 <= n )
    		{
	if (m==3)
		{difplot.lm(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE)
		}
	else 
		{difplot.lm2D(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,main=disti.sort[t1,1])
		}
      
      	cat(paste("outlier #",t1,": ",disti.sort[t1,1],"     procrustes dist. to mean: ",disti.sort[t1,2],"\n",sep=""))
      	
	if (disti.sort[t1,1] %in% outlier)
		{ answer<-substr(readline(" already added to outlierlist! remove fromlist (y/N/s)?"), 1L,1L)
		while (!(answer %in% c("y","Y","s","S","N","n")))
			{answer<-substr(readline(" already added to outlierlist! remove fromlist (y/N/s)?"), 1L,1L)
			}
		### keep in outlier list		
		if (answer == "n" || answer == "N")
         		{
          		lo<-lo+1
          		}
		### switching landmark routine starting ###	 	
		else if (answer == "s" || answer == "S")  
			{loop0="y"
			while(loop0=="y" || loop0 =="Y")
				{answer1a<-as.integer(readline("select  two Landmarks to switch position\ninsert first: "))
				while (is.na(answer1a)|| answer1a > k || answer1a < 1)
					{answer1a<-as.integer(readline(paste("please enter integer <",k,"!\ninsert first: ")))
					}
				answer1b<-as.integer(readline("insert second: "))
				while (is.na(answer1b) || answer1b > k || answer1b < 1)
					{answer1b<-as.integer(readline(paste("please enter integer <",k,"!\ninsert second: ")))
					}
		
		### switch rows of selected landmarks ##
			raw[c(answer1a,answer1b),,disti.sort[t1,1]]<-raw[c(answer1b,answer1a),,disti.sort[t1,1]]
			A$rotated[c(answer1a,answer1b),,disti.sort[t1,1]]<-A$rotated[c(answer1b,answer1a),,disti.sort[t1,1]]
			rho.new<-angle.calc(A$rotated[,,disti.sort[t1,1]],A$mshape)$rho
			
			if (m==3)
				{rgl.clear()
				difplot.lm(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE)
				}
			else 
				{difplot.lm2D(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,main=disti.sort[t1,1])
				}
			cat(paste("new distance to mean:",rho.new,"\n"))
			loop0<-substr(readline("switch more (y/N)? "),1L,1L)
				while(loop0 != "y" && loop0 != "Y" && loop0 != "n" && loop0 != "N" )
					{loop0<-substr(readline("yes or no? "),1L,1L)
					}
			}
			answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
       				while (!(answer %in% c("y","Y","N","n")))
				{answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
				}
       	
		
			}
		### remove from outlier list ###
		else if (answer == "y" || answer == "Y")
         			{outlier<-outlier[-which(outlier==disti.sort[t1,1])]
				cat(paste("outlier #",t1,": ",disti.sort[t1,1],"removed from outlierlist.\n"))
          			lo<-lo+1
          			}	
		}
      	
	else
	{	
	answer <- substr(readline(" add to outlierlist (y/N/s)?\ny=yes,n=no,s=switch landmarks: "), 1L,1L)
       	while (!(answer %in% c("y","Y","s","S","N","n")))
      		{answer <- substr(readline(" add to outlierlist (y/N/s)?\ny=yes,n=no,s=switch landmarks: "), 1L,1L)
       		}
	### add to outlier list		
		if (answer == "y" || answer == "Y")
         		{outlier[lo]<-disti.sort[t1,1]
          		lo<-lo+1
          		}
	### switching landmark routine starting ###	 	
		else if (answer == "s" || answer == "S")  
			{loop0="y"
			while(loop0=="y" || loop0 =="Y")
				{answer1a<-as.integer(readline("select  two Landmarks to switch position\ninsert first: "))
				while (is.na(answer1a)|| answer1a > k || answer1a < 1)
					{answer1a<-as.integer(readline(paste("please enter integer <",k,"!\ninsert first: ")))
					}
				answer1b<-as.integer(readline("insert second: "))
				while (is.na(answer1b) || answer1b > k || answer1b < 1)
					{answer1b<-as.integer(readline(paste("please enter integer <",k,"!\ninsert second: ")))
					}
		
		### switch rows of selected landmarks ##
			raw[c(answer1a,answer1b),,disti.sort[t1,1]]<-raw[c(answer1b,answer1a),,disti.sort[t1,1]]
			A$rotated[c(answer1a,answer1b),,disti.sort[t1,1]]<-A$rotated[c(answer1b,answer1a),,disti.sort[t1,1]]
			rho.new<-angle.calc(A$rotated[,,disti.sort[t1,1]],A$mshape)$rho
			if (m==3)
				{rgl.clear()
				difplot.lm(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE)
				}
			else 
				{difplot.lm2D(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,main=disti.sort[t1,1])
				}
			cat(paste("new distance to mean:",rho.new,"\n"))
			loop0<-substr(readline("switch more (y/N)? "),1L,1L)
				while(loop0 != "y" && loop0 != "Y" && loop0 != "n" && loop0 != "N" )
					{loop0<-substr(readline("yes or no? "),1L,1L)
					}
			}
			answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
       				while (!(answer %in% c("y","Y","N","n")))
				{answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
				}
			if (answer == "y" || answer == "Y")
         			{outlier[lo]<-disti.sort[t1,1]
          			lo<-lo+1
          			}	
			}
		    
	}	
      ### go on or halt ###
      	answer1 <- substr(readline(" next/previous/stop (n/p/s)?  "), 1L,1L)
	 	if (!(answer1 %in% c("P","p","S","s")))
		
			{rgl.clear()
           		rgl.bg(color = "white")
            		t1<-t1+1          
         		}
		else if (answer1 %in% c("P","p"))
			{if(m==3)
				{rgl.clear()
           			rgl.bg(color = "white")
				}
			if (t1 !=1)            		
				{t1<-t1-1
				}
			else 
				{cat("already at top of the line!\n")
				}        
         		}
			
    
		else if (answer1 == "s" || answer1 == "S")
			{break}
      
     
    	}
	### remove outlier from array
	if (! is.null(outlier))
    		{raw<-raw[,,-outlier]
		}
	
    
    
    
    	return(list(data.cleaned=raw,outlier=outlier,dist.sort=disti.sort))
    
}
