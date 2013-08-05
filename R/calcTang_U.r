.calcTang_U <- function(datamatrix,SMvector,outlines,deselect=FALSE)
{     	
	if (is.list(outlines)==FALSE)
      		{
		outlines <- list(outlines)
		}
      
      	dims <- dim(datamatrix)[2]
      	k <- dim(datamatrix)[1]
      	if (deselect==TRUE)
      		{
		SMvector <- c(1:k)[-SMvector]
		}
      	m <- length(SMvector)
      	tanvec <- matrix(0,k,dims)
      	U <- matrix(0,dims*k,m)
      	Gamma0 <- c(datamatrix)
      
      
      
      	for ( j in 1:length(outlines))
        	{
		lt <- length(outlines[[j]])
        	temp <- outlines[[j]]

### procedure for open curves ####        	
		if (outlines[[j]][1]!= outlines[[j]][lt])
        		{
			for (i in 1:lt)
          
          			{
				if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt)
           				{
					tanvec[temp[i],] <- (datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
           				}
          
           
          			else if (temp[i]%in%SMvector==TRUE && i==1)
           				{
					tanvec[temp[i],] <- (datamatrix[temp[i],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i],]-datamatrix[temp[i+1],])^2))
           				}
          			else if (temp[i]%in%SMvector==TRUE && i==lt)
           				{
					tanvec[temp[i],] <- (datamatrix[temp[i-1],]-datamatrix[temp[i],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i],])^2))
          				}
          #else {tanvec[i,] <- c(rep(0,dims))}
          			} 
			}
        
### procedure for closed curves ####
        	else if (outlines[[j]][1]== outlines[[j]][lt]) 
          		{
			for (i in 1:(lt-1))
                    		{
				if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt)
           				{
					tanvec[temp[i],] <- (datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
           				}
          
           
         			 else if (temp[i]%in%SMvector==TRUE && i==1)
           			{
				tanvec[temp[i],] <- (datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])^2))
           			}
          
          #else {tanvec[i,] <- c(rep(0,dims))}
          			} 
			}
        	}
    
    	SMsort <- sort(SMvector)
    	for (i in 1:m)
      		{
		U[SMsort[i],i] <- tanvec[SMsort[i],1]
        	U[k+SMsort[i],i] <- tanvec[SMsort[i],2] 
        	if (dims==3)
        		{
        		U[2*k+SMsort[i],i] <- tanvec[SMsort[i],3]
        		}
      		}
    
    
    
    return(list(tanvec=tanvec,SMvector=SMvector,U=U,Gamma0=Gamma0))             
    
}     
