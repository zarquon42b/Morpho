procSym <- function(dataarray, scale=TRUE, reflect=TRUE, CSinit=TRUE,  orp=TRUE, tol=1e-05, pairedLM=NULL, sizeshape=FALSE, use.lm=NULL, center.part=FALSE, distfun=c("riemann", "angle"), SMvector=NULL, outlines=NULL, deselect=FALSE, recursive=TRUE,iterations=0, initproc=FALSE)
{	t0 <- Sys.time()     
	A <- dataarray
      	k <- dim(A)[1]
      	m <- dim(A)[2]     
      	n <- dim(A)[3]
        Mir <- rep(1, m)
        Mir[1] <- -1
      	Mir <- diag(Mir)
      	dataslide <- NULL
      	CS <- NULL
        if (substr(distfun[1], 1L, 1L) == "r")
            distfun <- kendalldist
        else
            distfun <- function(x,y)
                {
                    rho <- angle.calc(x,y)
                    return(rho)
                }
        
        if (is.null(SMvector))   
            { 
                CS <- apply(A, 3, cSize)
            }
        
        if (!is.null(SMvector)) # includes sliding of Semilandmarks
            {
                if (is.null(outlines))
                    stop("please specify outlines")
                dataslide <- Semislide(A, SMvector=SMvector,outlines=outlines,tol=tol,deselect=deselect,recursive=recursive,iterations=iterations,pairedLM=pairedLM,initproc=initproc)
                A <- dataslide
                for (i in 1:n)
                    CS <- apply(A,3,cSize)
                if (CSinit==TRUE)
                    { 
                        for (i in 1:n)
                            A[,,i] <- A[,,i]/CS[i]
                    }
            }
###### create mirrored configs ######
        if (!is.null(pairedLM))
            {
                Amir <- A
                for (i in 1:n)
                    {
                        Amir[,,i] <- A[,,i]%*%Mir
                        Amir[c(pairedLM),,i] <- Amir[c(pairedLM[,2:1]),,i]
                    }
                Aall <- bindArr(A,Amir,along=3)
            }
        else
            Aall <- A
                
###### proc fit of all configs ######
      	cat("performing Procrustes Fit ")
        
        if (!is.null(use.lm)) ### only use subset for rotation and scale
            {
                proc <- ProcGPA(Aall[use.lm,,],scale=scale,CSinit=CSinit,reflection=reflect)
                tmp <- Aall
                for (i in 1:dim(Aall)[3])
                    {
                        tmp[,,i] <- rotonmat(Aall[,,i],Aall[use.lm,,i],proc$rotated[,,i],scale=TRUE, reflection=reflect)
                        if (center.part)
                            tmp[,,i] <- apply(tmp[,,i],2,scale,scale=F) ## center shapes
                        else
                            orp <- FALSE
                    }
                proc$rotated <- tmp
                proc$mshape <- apply(tmp,1:2,mean) ##calc new meanshape
            }
        else
            proc <- ProcGPA(Aall,scale=scale,CSinit=CSinit, reflection=reflect)
        
        procrot <- proc$rotated
        dimna <- dimnames(dataarray)
        if (!is.null(pairedLM))
            dimna[[3]] <- c(dimna[[3]],dimna[[3]])
        dimnames(proc$rotated) <- dimna
       	
 	meanshape <- proc$mshape
        rho <- NULL
        
        for (i in 1:n)
            rho[i] <- distfun(proc$rotated[,,i],proc$mshape)
        
        rmsrho <- sqrt(mean(rho^2))
        orpdata <- 0

###### project into tangent space ######
###test###        
                                        #meanall <- apply(proc$rotated[,,1:n],c(1,2),mean)
        if (orp==TRUE && CSinit==TRUE)
            procrot <- orp(proc$rotated, mshape=proc$mshape)
        
        orpdata <- procrot
        dimnames(orpdata) <- dimna
        
###### calculate Symmetric means ######
        if (!is.null(pairedLM)) 
            {
### generate symmetrized mean for each individual between original and mirrored configuration ###      		
                Symarray <- A
                for (i in 1:n)
                    Symarray[,,i] <- (procrot[,,i]+procrot[,,n+i])/2
### generate deviation between each individual and its specific symmetrized mean ###    		
                Asymm <- A 
                for (i in 1:n)
                    Asymm[,,i] <- (procrot[,,i]-Symarray[,,i])
                
                dimnames(Asymm) <- dimnames(dataarray)
            }
	else 
            Symarray <- procrot

        Symtan <- Symarray
        for (i in 1:n)
            Symtan[,,i] <- Symarray[,,i]-meanshape
        
      	tan <- matrix(NA,n,m*k)
      	for(i in 1:n)
            tan[i,] <- c(Symtan[,,i])
        
	if (sizeshape==TRUE)
            { 
                CSlog <- log(CS)-mean(log(CS))
                tan <- cbind(CSlog,tan)
            }
      	dimnames(Symarray) <- dimnames(dataarray)
        
###### PCA Sym Component ###### 
        princ <- try(prcomp(tan),silent=TRUE)
        if (class(princ) == "try-error")
            princ <- eigenPCA(tan)

        values <- 0
      	eigv <- princ$sdev^2
	
       	for (i in 1:length(eigv))
            {
                if (eigv[i] > 1e-14)
                    values[i] <- eigv[i]
            }
	lv <- length(values)
	PCs <- princ$rotation[,1:lv]
 	PCscore_sym <- as.matrix(princ$x[,1:lv])
        rownames(PCscore_sym) <- dimnames(dataarray)[[3]]
	rownames(tan) <- rownames(PCscore_sym)

###### create a neat variance table for Sym ###### 
        if (length(values)==1)
            SymVar <- values
        else
            {
                SymVar <- matrix(NA,length(values),3)
                SymVar[,1] <- values
                
                for (i in 1:length(values))
                    SymVar[i,2] <- (values[i]/sum(values))*100
                SymVar[1,3] <- SymVar[1,2]
                
                for (i in 2:length(values))
                    SymVar[i,3] <- SymVar[i,2]+ SymVar[i-1,3]
                
                colnames(SymVar) <- c("eigenvalues","% Variance","Cumulative %")
            }
        
        
###### PCA Asym Component ###### 
      	asvalues <- 0
      	PCs_Asym <- 0
        if (!is.null(pairedLM)) 
            {
                asymmean <- apply(Asymm,c(1,2),mean)
                asymtan <- matrix(NA,n,m*k)
                for(i in 1:n)
                    asymtan[i,] <- c(Asymm[,,i]-asymmean)
                
                pcasym <- try(prcomp(asymtan),silent=TRUE)
                if (class(pcasym) == "try-error")
                    pcasym <- eigenPCA(asymtan)
                
                asvalues <- 0
                eigva <- pcasym$sdev^2
                for (i in 1:length(eigva))
                    {
                        if (eigva[i] > 1e-14)
                            asvalues[i] <- eigva[i]
                    }
                lva <- length(asvalues)
                PCs_Asym <- pcasym$rotation[,1:lva]
                PCscore_asym <- as.matrix(pcasym$x[,1:lva])
                rownames(PCscore_asym) <- dimnames(dataarray)[[3]]
                rownames(asymtan) <- rownames(PCscore_sym)
                
###### create a neat variance table for Asym ######
                if (length(asvalues)==1)
                    AsymVar <- asvalues
                else
                    {
                        AsymVar <- matrix(NA,length(asvalues),3)
                        AsymVar[,1] <- asvalues
                        
                        for (i in 1:length(asvalues))
                            AsymVar[i,2] <- (asvalues[i]/sum(asvalues))*100
                        
                        AsymVar[1,3] <- AsymVar[1,2]
                        for (i in 2:length(asvalues))
                            AsymVar[i,3] <- AsymVar[i,2]+ AsymVar[i-1,3]
                        
                        colnames(AsymVar) <- c("eigenvalues","% Variance","Cumulative %")
                    }
            } 
###### output ######
        t1 <- Sys.time()
	cat(paste("Operation completed in",t1-t0,"secs\n"))
        
	
	if (!is.null(pairedLM))
            {
                out <- (list(
                    size=CS,rotated=proc$rotated[,,1:n],
                    rotmir=proc$rotated[,,(n+1):(2*n)],
                    Sym=Symarray,Asym=Asymm,asymmean=asymmean,
                    mshape=(meanshape+asymmean), symmean=meanshape,
                    Symtan=tan,Asymtan=asymtan,PCsym=PCs,PCscore_sym=PCscore_sym,
                    eigensym=values,SymVar=SymVar,PCasym=PCs_Asym,
                    PCscore_asym=PCscore_asym,eigenasym=asvalues,
                    AsymVar=AsymVar,orpdata=orpdata[,,1:n],
                    orpmir=orpdata[,,(n+1):(2*n)],
                    rmsrho=rmsrho,rho=rho,dataslide= dataslide,
                    pairedLM=pairedLM
                    ))
             class(out) <- "symproc"
             return(out)
         }
        
      	else  {
            out <- (list(
            size=CS,rotated=proc$rotated,mshape=meanshape,tan=tan,PCs=PCs,
            PCscores=PCscore_sym,eigenvalues=values,Variance=SymVar,
            orpdata=orpdata[,,1:n] ,rmsrho=proc$rmsrho,rho=rho,
            dataslide= dataslide
            ))
               class(out) <- "nosymproc"
               return(out)
           }
    }
