calcTang_U_s<-function(datamatrix,normalmatrix=NULL,SMvector,outlines=NULL,surface=NULL,deselect=FALSE)
{     	
	
		
	
      
      	dims<-dim(datamatrix)[2]
      	k<-dim(datamatrix)[1]
    	if (deselect==TRUE)
      		{
		SMvector<-c(1:k)[-SMvector]
		}
      	m<-length(SMvector)
	if (is.null(surface))
		{tanvec<-matrix(0,k,dims)
      		U<-matrix(0,dims*k,m)
		}
	else
		{tanvec<-matrix(0,k,dims*2)
      		U<-matrix(0,dims*k,m*2)
		}
      	Gamma0<-c(datamatrix)
      
      
	if (is.null(outlines) == FALSE)      
    
		{  			
			if (is.list(outlines)==FALSE)
      				{
				outlines<-list(outlines)
				}
			for ( j in 1:length(outlines))
        		{
			lt<-length(outlines[[j]])
        		temp<-outlines[[j]]

### procedure for open curves ####        	
		if (outlines[[j]][1]!= outlines[[j]][lt])
        		{
			for (i in 1:lt)
          
          			{
				if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt)
           				{
					tanvec[temp[i],1:3]<-(datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
           				}
          
           
          			else if (temp[i]%in%SMvector==TRUE && i==1)
           				{
					tanvec[temp[i],1:3]<-(datamatrix[temp[i],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i],]-datamatrix[temp[i+1],])^2))
           				}
          			else if (temp[i]%in%SMvector==TRUE && i==lt)
           				{
					tanvec[temp[i],1:3]<-(datamatrix[temp[i-1],]-datamatrix[temp[i],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i],])^2))
          				}
          #else {tanvec[i,]<-c(rep(0,dims))}
          			} 
			}
        
### procedure for closed curves ####
        	else if (outlines[[j]][1]== outlines[[j]][lt]) 
          		{
			for (i in 1:(lt-1))
                    		{
				if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt)
           				{
					tanvec[temp[i],1:3]<-(datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
           				}
          
           
         			 else if (temp[i]%in%SMvector==TRUE && i==1)
           			{
				tanvec[temp[i],1:3]<-(datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])^2))
           			}
          
          #else {tanvec[i,]<-c(rep(0,dims))}
          			} 
			}
        	}
	}
    
### procedure for surfaces ###
	if (is.null (surface) ==F)
	 	    	{	
			lt<-length(surface)
       			temp<-surface
			for (i in 1:lt)
                    		{
					tanp<-tanplan(normalmatrix[temp[i],])
					tanvec[temp[i],]<-c(tanp$y,tanp$z)				
				}
			}
#		pointset<-datamatrix[temp,]
#		write.obj(cbind("v",pointset),filename="temp")
#		write(paste("<!DOCTYPE FilterScript>\n","<FilterScript>\n"," <filter name=\"Compute normals for point sets\">\n", "<Param type=\"RichInt\" value=\"10\" name=\"K\"/>\n","</filter>\n","</FilterScript>",sep=""),file="norm.mlx")
#		}
#		command<-"meshlabserver -i temp.obj -o temp1.obj -s norm.mlx -om vn"
#		system(command)


#### end surfaces ####
	    	
	SMsort<-sort(SMvector)
	if (!is.null(surface))    	
		{		
			for (i in 1:m)
      			{
			U[SMsort[i],i]<-tanvec[SMsort[i],1]
        		U[k+SMsort[i],i]<-tanvec[SMsort[i],2] 
        		U[2*k+SMsort[i],i]<-tanvec[SMsort[i],3]
  			U[SMsort[i],(i+m)]<-tanvec[SMsort[i],4]
			U[k+SMsort[i],(i+m)]<-tanvec[SMsort[i],5]
			U[2*k+SMsort[i],(i+m)]<-tanvec[SMsort[i],6]
			}
		}
	
	else 	
		{
			for (i in 1:m)
      			{
			U[SMsort[i],i]<-tanvec[SMsort[i],1]
        		U[k+SMsort[i],i]<-tanvec[SMsort[i],2] 
        		U[2*k+SMsort[i],i]<-tanvec[SMsort[i],3]
			}
		}
			
    
    
    
    return(list(tanvec=tanvec,SMvector=SMvector,U=U,Gamma0=Gamma0))             
    
}     
