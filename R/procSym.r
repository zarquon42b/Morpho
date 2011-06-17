procSym<-function(dataarray,pairedLM=0,SMvector=0,outlines=0,orp=TRUE,tol=1e-05,CSinit=TRUE,deselect=FALSE,recursive=TRUE,iterations=0,scale=TRUE,reflect=FALSE,sizeshape=FALSE,initproc=FALSE,ignore.ref=FALSE)
{     
	A<-dataarray
      	k<-dim(A)[1]
      	m<-dim(A)[2]     
      	n<-dim(A)[3]
      	Mir<-diag(c(-1,1,1))
      	dataslide<-NULL
      	CS<-NULL
      
      	if (SMvector[1]==0) 
		{ 
      		CS<-apply(A,3,c.size)
		
		if (CSinit==TRUE)
          		{ 
			for (i in 1:n)
            			{A[,,i]<-A[,,i]/CS[i]}
          		}
		}
      		
      
      
      
      if (SMvector[1]!=0)           # includes sliding of Semilandmarks
      		{ 	
		dataslide<-Semislide(A, SMvector=SMvector,outlines=outlines,tol=tol,deselect=deselect,recursive=recursive,iterations=iterations,pcaoutput=FALSE,pairedLM=pairedLM,initproc=initproc)
        	A<-dataslide
        
        
        	for (i in 1:n)
        		CS<-apply(A,3,c.size)
        	if (CSinit==TRUE)
          		{ 
			for (i in 1:n)
            			{A[,,i]<-A[,,i]/CS[i]}
          		}
 	       }
###### create mirrored configs ######
        if (pairedLM[1]!=0)
        {
            Amir<-A
            for (i in 1:n)
              {Amir[,,i]<-A[,,i]%*%Mir
            Amir[c(pairedLM),,i]<-Amir[c(pairedLM[,2:1]),,i]}
            Aall<-abind(A,Amir)
        }
      
      else {Aall<-A}

###### proc fit of all configs ###### 
        proc<-procGPA(Aall,pcaoutput=FALSE,scale=scale,reflect=reflect)
        procrot<-proc$rotated
	
	dimna<-dimnames(dataarray)
	if (pairedLM[1]!=0)
		{			
			dimna[[m]]<-c(dimna[[m]],dimna[[m]])
		}
			
        	dimnames(proc$rotated)<-dimna
       	
 	meanshape<-proc$mshape


        
	rho<-proc$rho
          if (reflect==TRUE && ignore.ref==TRUE)
          {
          for (i in 1:n)
            {rho[i]<-riemdist(proc$rotated[,,i],proc$mshape)}
          }
         
        orpdata<-0

###### project into tangent space ######
	###test###        
		#meanall<-apply(proc$rotated[,,1:n],c(1,2),mean)
        if (orp==TRUE && CSinit=TRUE)
        	{
		procrot<-orp(proc$rotated)
        	}
		orpdata<-procrot
        	dimnames(orpdata)<-dimna
		
      
###### calculate Symmetric means ######
	if (pairedLM[1]!=0) 
      		{
 	### generate symmetrized mean for each individual between original and mirrored configuration ###      		
		Symarray<-A
      		for (i in 1:n)
      			{
			Symarray[,,i]<-(procrot[,,i]+procrot[,,n+i])/2
      			}
  	### generate deviation between each individual and its specific symmetrized mean ###    		
		Asymm<-A 
      		for (i in 1:n)
      			{
			Asymm[,,i]<-(procrot[,,i]-Symarray[,,i])
      			}
      		dimnames(Asymm)<- dimnames(dataarray)
      		}
	else 
      		{
      		Symarray<-procrot
      		}
      
	Symtan<-Symarray
      
	for (i in 1:n)
      		{Symtan[,,i]<-Symarray[,,i]-meanshape}
      
      	tan<-matrix(NA,n,m*k)
      	for(i in 1:n)
      		{tan[i,]<-c(Symtan[,,i])}
        
	if (sizeshape==TRUE)
          	{ 
		CSlog<-log(CS)-mean(log(CS))
            	tan<-cbind(CSlog,tan)
          	}
      	dimnames(Symarray)<-dimnames(dataarray)
      
###### PCA Sym Component ###### 
       		princ<-prcomp(tan)
	values<-0
      	eigv<-princ$sdev^2
	
       	for (i in 1:length(eigv))
       		{
		if (eigv[i] > 1e-14)
        		{
			values[i]<-eigv[i]
         		}
        	}
	lv<-length(values)
	PCs<-princ$rotation[,1:lv]
 	PCscore_sym<-princ$x[,1:lv]

###### create a neat variance table for Sym ###### 
        if (length(values)==1)
          	{SymVar<-values}
        else
        	{
          	SymVar<-matrix(NA,length(values),3)
          	SymVar[,1]<-values
        
          	for (i in 1:length(values))
            		{
              		SymVar[i,2]<-(values[i]/sum(values))*100
            		}
          	SymVar[1,3]<- SymVar[1,2]
          	for (i in 2:length(values))
           		{         
             		SymVar[i,3]<-SymVar[i,2]+ SymVar[i-1,3]
            		}
          	colnames(SymVar)<-c("eigenvalues","% Variance","Cumulative %")
        	}
      
      
###### PCA Asym Component ###### 
      	asvalues<-0
      	PCs_Asym<-0
      	if (pairedLM[1]!=0) 
      		{
      		asymtan<-matrix(NA,n,m*k)
      		for(i in 1:n)
      			{ 
       		 	asymmean<-apply(Asymm,c(1,2),mean)
        		asymtan[i,]<-c(Asymm[,,i]-asymmean)
			}
       
      		pcasym<-prcomp(asymtan)
       		asvalues<-0
       		eigva<-pcasym$sdev^2
		for (i in 1:length(eigv))
       		{
			if (eigva[i] > 1e-14)
        			{
				asvalues[i]<-eigva[i]
         			}
        		}
        	lva<-length(asvalues)
		PCs_Asym<-pcasym$rotation[,1:lva]
         	PCscore_asym<-pcasym$x[,1:lva]


###### create a neat variance table for Asym ######
        	if (length(asvalues)==1)
          		{
			AsymVar<-asvalues
			}
        	else
        		{
         		AsymVar<-matrix(NA,length(asvalues),3)
          		AsymVar[,1]<-asvalues
        
          		for (i in 1:length(asvalues))
            			{
              			AsymVar[i,2]<-(asvalues[i]/sum(asvalues))*100
            			}
          		AsymVar[1,3]<- AsymVar[1,2]
          		for (i in 2:length(asvalues))
            			{         
              			AsymVar[i,3]<-AsymVar[i,2]+ AsymVar[i-1,3]
            			}
          		colnames(AsymVar)<-c("eigenvalues","% Variance","Cumulative %")
        		}
		} 
###### output ######
	
	if (pairedLM[1]!=0)
      	{return(list(size=CS,rotated=proc$rotated[,,1:n],rotmir=proc$rotated[,,(n+1):(2*n)],Sym=Symarray,Asym=Asymm,asymmean=asymmean,mshape=(meanshape+asymmean),
	symmean=meanshape,Symtan=tan,Asymtan=asymtan,PCsym=PCs,PCscore_sym=PCscore_sym,eigensym=values,SymVar=SymVar,PCasym=PCs_Asym,PCscore_asym=PCscore_asym,eigenasym=asvalues,AsymVar=AsymVar,orpdata=orpdata[,,1:n],orpmir=orpdata[,,(n+1):(2*n)],rmsrho=proc$rmsrho,rho=rho,dataslide= dataslide))
      }
      
      	else  {return(list(size=CS,rotated=proc$rotated,mshape=meanshape,tan=tan,PCs=PCs,PCscores=PCscore_sym,eigenvalues=values,Variance=SymVar,orpdata=orpdata[,,1:n] ,rmsrho=proc$rmsrho,rho=rho,dataslide= dataslide))
      }
}
