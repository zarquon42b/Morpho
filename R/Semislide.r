Semislide <- function(dataframe,SMvector,outlines,tol=1e-05,deselect=FALSE,recursive=TRUE,iterations=0,initproc=FALSE,pairedLM=NULL,bending=TRUE,stepsize=1)
{
    n <- dim(dataframe)[3]
    k <- dim(dataframe)[1]
    m <- dim(dataframe)[2]
    p1 <- 10^12
    if (iterations == 0)
        iterations <- 1e10
    
    ini <- rotonto(dataframe[,,1],dataframe[,,2],reflection=T)
    mshape <- (ini$X+ini$Y)/2
    
    if(initproc==TRUE) { # perform proc fit before sliding
        procini <- ProcGPA(dataframe,scale=TRUE)
        mshape <- procini$mshape
    }
    dataslide <- dataframe
    
    if (!is.null(pairedLM)) {# create symmetric mean to get rid of assymetry along outline after first relaxation
        Mir <- diag(c(-1,1,1))
        A <- mshape
        Amir <- mshape%*%Mir
        Amir[c(pairedLM),] <- Amir[c(pairedLM[,2:1]),]
        symproc <- rotonto(A,Amir)
        mshape <- (symproc$X+symproc$Y)/2
    }
    count <- 1
    while (p1>tol && count <= iterations) {
        dataslide_old <- dataslide
        mshape_old <- mshape
        cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated  
        
        if (recursive==TRUE)      # slided Semilandmarks are used in next iteration step
            dataframe <- dataslide
        
        L <- CreateL(mshape)
        
        for (j in 1:n) {
            tmp <- dataframe[,,j]
            if (!bending) {
                rot <- rotonto(mshape,tmp,scale=TRUE,reflection=FALSE)
                tmp <- rot$yrot
            }
            U <- .calcTang_U(tmp,SMvector=SMvector,outlines=outlines,deselect=deselect)
            if (bending) {
                dataslide[,,j] <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m,stepsize=stepsize)
            } else {
                tmpslide <- calcProcDGamma(U$U,U$Gamma0,mshape,dims=m,stepsize=stepsize)
                dataslide[,,j] <- rotreverse(tmpslide,rot)
            }
        }
        proc <- ProcGPA(dataslide,scale=TRUE)
        mshape <- proc$mshape
        p1_old <- p1   
        p1 <- sum(diag(crossprod((mshape_old/cSize(mshape_old))-(mshape/cSize(mshape)))))
                       
        ## check for increasing convergence criterion ###		
        if (p1 > p1_old) {
            dataslide <- dataslide_old
            cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
            p1 <- 0
        } else {
            cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
            count <- count+1 
        }          		
    }
    return(dataslide)
}
