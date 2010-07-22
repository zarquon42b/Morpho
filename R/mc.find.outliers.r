mc.find.outliers<-function(A,color=4,lwd=1,lcol=2)
{   	
	raw<-A
	n<-dim(A)[3]
	k<-dim(A)[1]
	A<-mc.procGPA(A)
	
	rho<-0
	for(i in 1:n)
		{rho[i]<-angle.calc(A$rotated[,,i],A$mshape)$rho
		}
	A$rho<-rho
	
	if (dim (A$rotated)[2]==2)
		{A0<-array(NA,dim=c(k,3,n))
		for (i in 1:n)
			
			{A0[,,i]<-cbind(A$rotated[,,i],0)
			}
		A$mshape<-cbind(A$mshape,0)
		A$rotated<-A0
		}
    	#rgl.clear()
    	#rgl.bg(color = "white")
   	disti<-data.frame(c(1:n),A$rho)
    	disti.sort<-disti[order(disti[,2],decreasing=T),]
    	colnames(disti.sort)[1]<-"# in array"
    	rownames(disti.sort)<-c(1:n)
    	outlier<-NULL
    	t1 <- 1
    	lo<-1
    	
    	while (t1 <= n )
    		{
		difplot.lm(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE)
      
      	cat(paste("outlier #",t1,": ",disti.sort[t1,1],"     procrustes dist. to mean: ",disti.sort[t1,2],"\n",sep=""))
      
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
				while (is.na(answer1a))
					{answer1a<-as.integer(readline("please enter integer!\ninsert first: "))
					}
				answer1b<-as.integer(readline("select  two Landmarks to switch position\nsecond: "))
				while (is.na(answer1b))
					{answer1b<-as.integer(readline("please enter integer!\ninsert second: "))
					}
		
		### switch rows of selected landmarks ##
			raw[c(answer1a,answer1b),,disti.sort[t1,1]]<-raw[c(answer1b,answer1a),,disti.sort[t1,1]]
			A$rotated[c(answer1a,answer1b),,disti.sort[t1,1]]<-A$rotated[c(answer1b,answer1a),,disti.sort[t1,1]]
			rho.new<-angle.calc(A$rotated[,,disti.sort[t1,1]],A$mshape)$rho
			rgl.clear()
			difplot.lm(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE)
			cat(paste("new distance to mean:",rho.new,"\n"))
			loop0<-substr(readline("switch more (y/N)? "),1L,1L)
				if(loop0 != "y" && loop0 != "Y" && loop0 != "n" && loop0 != "N" )
					{loop0<-substr(readline("yes or no? "),1L,1L)
					}
			}
			answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
       	
			if (answer == "y" || answer == "Y")
         			{outlier[lo]<-disti.sort[t1,1]
          			lo<-lo+1
          			}	
			}
    
		
      ### go on or halt ###
      	answer1 <- substr(readline(" wanna go on? (y/N)?  "), 1L,1L)
	 	if (!(answer1 %in% c("N","Y","n","y")))
		
		{rgl.clear()
           	rgl.bg(color = "white")
            t1<-t1+1          
         	}
    
		else if (answer1 == "n" || answer1 == "N")
			{break}
      
     
    	}
	### remove outlier from array
	if (! is.null(outlier))
    		{raw<-raw[,,-outlier]
		}
	
    
    
    
    	return(list(data.cleaned=raw,outlier=outlier,dist.sort=disti.sort))
    
}
