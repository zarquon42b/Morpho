Semislide<-function(dataframe,SMvector,outlines,tol=1e-05,deselect=FALSE,recursive=TRUE,iterations=0,pcaoutput=TRUE,initproc=FALSE,pairedLM=0)

{     n<-dim(dataframe)[3]
      k<-dim(dataframe)[1]
      m<-dim(dataframe)[2]
      
      p1<-10^12
      
      ini<-procOPA(dataframe[,,1],dataframe[,,2],reflect=T)
	mshape<-(ini$Ahat+ini$Bhat)/2
      
      if(initproc==TRUE) # perform proc fit before sliding
      {procini<-procGPA(dataframe,pcaoutput=F,scale=TRUE,reflect=TRUE,distances=FALSE)
        mshape<-procini$mshape
        
      }
      dataslide<-dataframe
      
      if (pairedLM[1]!=0)# create symmetric mean to get rid of assymetry along outline after first relaxation
      {
      Mir<-diag(c(-1,1,1))
      A<-mshape
      Amir<-mshape%*%Mir
      Amir[c(pairedLM),]<-Amir[c(pairedLM[,2:1]),]
      symproc<-procOPA(A,Amir)
      mshape<-(A+Amir)/2
      }
      
      if (iterations!=0) # calculation for a defined number of iterations
                   
      		{count<-1
        	for (i in 1:iterations)
        		{dataslide_old<-dataslide
			mshape_old<-mshape
          		cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated  
        
        		if (recursive==TRUE)      # slided Semilandmarks are used in next iteration step
          			{ dataframe<-dataslide
				}
            		if (m==3)
            			{L<-CreateL(mshape)
				}
            		else 
            			{L<-CreateL2D(mshape)
				} 
          		for (j in 1:n)
		  		{U<-calcTang_U(dataframe[,,j],SMvector=SMvector,outlines=outlines,deselect=deselect)
            			dataslide[,,j]<-calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m)$Gamatrix
          			}
          		proc<-procGPA(dataslide,scale=TRUE,reflect=TRUE,pcaoutput=FALSE,distances=FALSE)
          		mshape<-proc$mshape
			p1_old<-p1   
			p1<-sum(diag(crossprod((mshape_old/c.size(mshape_old))-(mshape/c.size(mshape)))))
         		 #p1<-sum(diag(crossprod(mshape_old-mshape)))/k
			if (p1 > p1_old)
				{
				dataslide<-dataslide_old
				cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
				p1<-0
				} 
			else
				{
	          		cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
         	 		count<-count+1 
				}          		
			     
        		}
      		}
      
      
      	if (iterations==0)   # calculation until convergence of meanshape
      		{count<-1
        	while (p1>tol)
         		{dataslide_old<-dataslide 
			mshape_old<-mshape
           
           		cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated
          		if (recursive==TRUE)    # slided Semilandmarks are used in next iteration step
          			{dataframe<-dataslide
				}
            		if (m==3)
            			{L<-CreateL(mshape)
				}
            		else 
            			{L<-CreateL2D(mshape)
				} 
            		for (j in 1:n)
          			{U<-calcTang_U(dataframe[,,j],SMvector=SMvector,outlines=outlines,deselect=deselect)
            			dataslide[,,j]<-calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m)$Gamatrix
          			}
          		proc<-procGPA(dataslide,scale=TRUE,reflect=TRUE,pcaoutput=pcaoutput,,distances=FALSE)
          		mshape<-proc$mshape
          		p1_old<-p1   
			p1<-sum(diag(crossprod((mshape_old/c.size(mshape_old))-(mshape/c.size(mshape)))))
         		 #p1<-sum(diag(crossprod(mshape_old-mshape)))/k
			if (p1 > p1_old)
				{
				dataslide<-dataslide_old
				cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
				p1<-0
				} 
			else
				{
	          		cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
         	 		count<-count+1 
				}    
      			}
		}

        
        return(dataslide)
}


