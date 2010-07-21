\name{mc.procSym}
\alias{mc.procSym}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{performs a  Procrustes Analysis, including options for symmetrical configurations and sliding of 2D and 3D -curves.
  
}
\description{
    \code{procSym} performs Procrustes superimposition including sliding of semilandmarks in 3d
}
\usage{
mc.procSym(dataarray, pairedLM = 0, SMvector = 0, outlines = 0, orp = TRUE, tol = 1e-05, CSinit = TRUE, deselect = FALSE, recursive = TRUE, iterations = 0, scale = TRUE, reflect = FALSE, sizeshape = FALSE, initproc = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataarray}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
%%     ~~Describe \code{dataarray} here~~
}
  \item{pairedLM}{A X x 2 matrix containing the paired LM. E.g. the left column contains the lefthand landmarks, 
                  while the right side contains the corresponding right hand landmarks.
%%     ~~Describe \code{pairedLM} here~~
}
  \item{SMvector}{A vector containing the landmarks on the curve(s) that are allowed to slide
%%     ~~Describe \code{SMvector} here~~
}
  \item{outlines}{A vector (or if threre are several curves) a list of vectors (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in the successive position on the curve - including the beginning and end points, that are not allowed to slide.
%%     ~~Describe \code{outlines} here~~
}
  \item{orp}{Logical: if TRUE, an orthogonal projection into shape space is performed.
%%     ~~Describe \code{orp} here~~
}
  \item{tol}{numeric: Threshold for convergence in the sliding process
%%     ~~Describe \code{tol} here~~
}
  \item{CSinit}{Logical: if TRUE, all configurations are initially scaled to Unit Centroid Size.
%%     ~~Describe \code{CSinit} here~~
}
  \item{deselect}{Logical: if TRUE, the SMvector is interpreted as those landmarks, that are not allowed to slide.
%%     ~~Describe \code{deselect} here~~
}
  \item{recursive}{Logical: if TRUE, during the iterations of the sliding process, the outcome of the previous iteration will be used. 
                    Otherwise the original configuration will be used in all iterations. 
%%     ~~Describe \code{recursive} here~~
}
  \item{iterations}{integer: select manually how many iterations will be performed during the sliding process 
                    (usefull, when there is very slow convergence). 
                              0 means iteration until convergence.
%%     ~~Describe \code{iterations} here~~
}
  \item{scale}{Logical: indicating if scaling is required
%%     ~~Describe \code{scale} here~~
}
  \item{sizeshape}{Logical: if TRUE, a log transformed variable of Centroid Size will be added to the shapedata as first variable before performing the PCA. 
%%     ~~Describe \code{reflect} here~~
}
  \item{sizeshape}{Logical: indicating if PCA in shape-size-space is requested
%%     ~~Describe \code{reflect} here~~
}
  \item{initproc}{Logical: indicating if the first Relaxation step is performed against the mean of an initial Procrustes superimposition. 
                  Symmetrical configurations will be relaxed against a perfectly symmetrical mean.
%%     ~~Describe \code{initproc} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{size }{a vector containing the Centroid Size of the configurations}
  
  \item{rotated }{k x m x n array of the rotated configurations}
  
  \item{Sym }{k x m x n array of the Symmetrical component - only available for the "Symmetry"-Option}
  
  \item{Asym }{k x m x n array of the Asymmetrical component - only available for the "Symmetry"-Option}
  \item{asymmean }{k x m matrix of mean asymmetric deviation from symmetric mean}
  \item{mshape }{sample meanshape}
  \item{symmean }{meanshape of symmetrized configurations}
  \item{tan }{if orp=TRUE: Residuals in tangentspace else, Procrustes residuals - only available without the "Symmetrie"-Option}
  \item{PCs }{Principal Components - if sizeshape=TRUE, the first variable of the PCs is size information (as log transformed Centroid Size)}
  \item{PCsym }{Principal Components of the Symmetrical Component}
  \item{PCasym }{Principal Components of the Asymmetrical Component}
  \item{PCscores }{PC scores}
  \item{PCscore_sym }{PC scores of the Symmetrical Component}
  \item{PCscore_asym }{PC scores of the Asymmetrical Component}
  \item{eigenvalues }{eigenvalues of the Covariance matrix}
  \item{eigensym }{eigenvalues of the "Symmetrical" Covariance matrix}
  \item{eigenasym }{eigenvalues of the "Asymmetrical" Covariance matrix}
  \item{Variance }{Table of the explained Variance by the PCs}
  \item{SymVar }{Table of the explained "Symmetrical" Variance by the PCs}
  \item{AsymVar }{Table of the explained "Asymmetrical" Variance by the PCs}
  \item{orpdata }{k x m x n array of the rotated configurations projected into tangent space}
  \item{rho }{vector of Riemannian distance from the mean}
}
\references{
\code{Dryden IL, and Mardia KV. 1998. Statistical shape analysis. Chichester.}
\code{Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape analysis of symmetric structures: quantifying variation among individuals and asymmetry. Evolution 56(10):1909-1920.}
\code{Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. «Semilandmarks in Three Dimensions,» in Modern Morphometrics in Physical Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer Academic/Plenum Publishers.}


}
\author{
Stefan Schlager
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(dataarray,pairedLM=0,SMvector=0,outlines=0,orp=TRUE,tol=1e-05,CSinit=TRUE,deselect=FALSE,recursive=TRUE,iterations=0,scale=TRUE,reflect=FALSE,sizeshape=FALSE,initproc=FALSE)
{     	t0<-Sys.time()
	A<-dataarray
      	k<-dim(A)[1]
      	m<-dim(A)[2]     
      	n<-dim(A)[3]
      	Mir<-diag(c(-1,1,1))
      	dataslide<-NULL
      	CS<-NULL
      
      	if (SMvector[1]==0) 
		{ 
      		CS<-apply(A,3,centroid.size)
		
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
        		CS<-apply(A,3,centroid.size)
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
              {Amir[,,i]<-A[,,i]\%*\%Mir
            Amir[c(pairedLM),,i]<-Amir[c(pairedLM[,2:1]),,i]}
            Aall<-abind(A,Amir)
        }
      
      else {Aall<-A}

###### proc fit of all configs ###### 
	cat("performing Procrustes Fit ")        
	proc<-mc.procGPA(Aall,scale=scale)
        procrot<-proc$rotated
	
	dimna<-dimnames(dataarray)
	if (pairedLM[1]!=0)
		{			
			dimna[[m]]<-c(dimna[[m]],dimna[[m]])
		}
			
        	dimnames(proc$rotated)<-dimna
       	
 	meanshape<-proc$mshape


        
	rho<-NULL
          {
          for (i in 1:n)
            {rho[i]<-angle.calc(proc$rotated[,,i],proc$mshape)$rho}
          }
         
        orpdata<-0

###### project into tangent space ######
	###test###        
		#meanall<-apply(proc$rotated[,,1:n],c(1,2),mean)
        if (orp==TRUE)
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
       	pcsym<-eigen(cov(tan))
       	values<-0
       	t1<-dim(tan)[2]
       	for (i in 1:length(pcsym$values))
       		{
		if (pcsym$values[i] > 1e-14)
        		{
			values[i]<-pcsym$values[i]
         		}
        	}
        
	PCs<-matrix(NA,t1,length(values))
        for (i in 1:length(values))
        	{     
		PCs[,i]<-pcsym$vectors[,i]
        	}
        PCs<-as.matrix(PCs)
	PCscore_sym<-tan\%*\%PCs

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
          	colnames(SymVar)<-c("eigenvalues","\% Variance","Cumulative \%")
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
       
      		pcasym<-eigen(cov(asymtan))
       		asvalues<-0
       		for (i in 1:length(pcasym$values))
       			{
			if (pcasym$values[i] > 1e-14)
        			{
				asvalues[i]<-pcasym$values[i]
         			}
        		}
        
		PCs_Asym<-matrix(NA,k*m,length(asvalues))
       		for (i in 1:length(asvalues))
        		{     
			PCs_Asym[,i]<-pcasym$vectors[,i]
        		}
        	PCs_Asym<-as.matrix(PCs_Asym)
        	PCscore_asym<-asymtan\%*\%PCs_Asym

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
          		colnames(AsymVar)<-c("eigenvalues","\% Variance","Cumulative \%")
        		}
		} 
###### output ######
	
	t1<-Sys.time()
	cat(paste("Operation completed in",t1-t0,"secs\n"))

	if (pairedLM[1]!=0)
      	{return(list(size=CS,rotated=proc$rotated[,,1:n],rotmir=proc$rotated[,,(n+1):(2*n)],Sym=Symarray,Asym=Asymm,asymmean=asymmean,mshape=(meanshape+asymmean),
	symmean=meanshape,Symtan=tan,Asymtan=asymtan,PCsym=PCs,PCscore_sym=PCscore_sym,eigensym=values,SymVar=SymVar,PCasym=PCs_Asym,PCscore_asym=PCscore_asym,eigenasym=asvalues,AsymVar=AsymVar,orpdata=orpdata[,,1:n],orpmir=orpdata[,,(n+1):(2*n)],rmsrho=proc$rmsrho,rho=rho,dataslide= dataslide))
      }
      
      	else  {return(list(size=CS,rotated=proc$rotated,mshape=meanshape,tan=tan,PCs=PCs,PCscores=PCscore_sym,eigenvalues=values,Variance=SymVar,orpdata=orpdata[,,1:n] ,rmsrho=proc$rmsrho,rho=rho,dataslide= dataslide))
      }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
