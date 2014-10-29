.calcTang_U_s <- function(datamatrix,normalmatrix=NULL,SMvector,outlines=NULL,surface=NULL,free=NULL,deselect=FALSE,weights=NULL)
{
    
    dims <- dim(datamatrix)[2]
    k <- dim(datamatrix)[1]
    
    if (is.null(weights)){
        weights <- 1
    } else
        weights <- c(weights,weights,weights)
    
    if (deselect==TRUE)
        SMvector <- c(1:k)[-SMvector]

    SMvector <- unique(SMvector)
    m <- length(SMvector)
    if ( !is.null(free)) {
        tanvec <- matrix(0,k,dims*3)
        U <- matrix(0,dims*k,m*3)
    } else if(!is.null(surface)) {
        tanvec <- matrix(0,k,dims*2)
        U <- matrix(0,dims*k,m*2)
    } else {
        tanvec <- matrix(0,k,dims)
        U <- matrix(0,dims*k,m)
    }
    Gamma0 <- c(datamatrix)
    
    if (is.null(outlines) == FALSE) {  			
        if (is.list(outlines)==FALSE) {
            outlines <- list(outlines)
        }
        for ( j in 1:length(outlines)) {
            lt <- length(outlines[[j]])
            temp <- outlines[[j]]
            
### procedure for open curves ####        	
            if (outlines[[j]][1]!= outlines[[j]][lt]) {
                for (i in 1:lt) {
                    if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt) {
                        tanvec[temp[i],1:3] <- (datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
                    } else if (temp[i]%in%SMvector==TRUE && i==1) {
                        tanvec[temp[i],1:3] <- (datamatrix[temp[i],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i],]-datamatrix[temp[i+1],])^2))
                    } else if (temp[i]%in%SMvector==TRUE && i==lt) {
                        tanvec[temp[i],1:3] <- (datamatrix[temp[i-1],]-datamatrix[temp[i],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i],])^2))
                    }
                } 
            }
### procedure for closed curves ####
            else if (outlines[[j]][1]== outlines[[j]][lt]) {
                for (i in 1:(lt-1)) {
                    if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt) {
                        tanvec[temp[i],1:3] <- (datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
                    } else if (temp[i]%in%SMvector==TRUE && i==1) {
                        tanvec[temp[i],1:3] <- (datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])^2))
                    }
                } 
            }
        }
    }
### procedure for surfaces ###
    if (is.null (surface) ==F) {	
        lt <- length(surface)
        temp <- surface
        for (i in 1:lt) {
            tanp <- tanplan(normalmatrix[temp[i],])
            tanvec[temp[i],1:6] <- c(tanp$y,tanp$z)				
        }
    }
#### end surfaces ####
    
### procedure for free sliding points ####
    
    if (!is.null(free)) {
        lf <- length(free)
        tmp <- free
        for (i in 1:lf)
            tanvec[tmp[i],] <- c(1,0,0,0,1,0,0,0,1)
    }
### end free sliding ##
    gc() 	
    SMsort <- sort(SMvector)
    if (!is.null(free)) {		
        for (i in 1:m) {
            U[SMsort[i],i] <- tanvec[SMsort[i],1]
            U[k+SMsort[i],i] <- tanvec[SMsort[i],2] 
            U[2*k+SMsort[i],i] <- tanvec[SMsort[i],3]
            U[SMsort[i],(i+m)] <- tanvec[SMsort[i],4]
            U[k+SMsort[i],(i+m)] <- tanvec[SMsort[i],5]
            U[2*k+SMsort[i],(i+m)] <- tanvec[SMsort[i],6]
            U[SMsort[i],(i+2*m)] <- tanvec[SMsort[i],7]
            U[k+SMsort[i],(i+2*m)] <- tanvec[SMsort[i],8]
            U[2*k+SMsort[i],(i+2*m)] <- tanvec[SMsort[i],9]
        }
    } else if (!is.null(surface)) {		
        for (i in 1:m) {
            U[SMsort[i],i] <- tanvec[SMsort[i],1]
            U[k+SMsort[i],i] <- tanvec[SMsort[i],2] 
            U[2*k+SMsort[i],i] <- tanvec[SMsort[i],3]
            U[SMsort[i],(i+m)] <- tanvec[SMsort[i],4]
            U[k+SMsort[i],(i+m)] <- tanvec[SMsort[i],5]
            U[2*k+SMsort[i],(i+m)] <- tanvec[SMsort[i],6]
        }
    } else {
        for (i in 1:m) {
            U[SMsort[i],i] <- tanvec[SMsort[i],1]
            U[k+SMsort[i],i] <- tanvec[SMsort[i],2] 
            U[2*k+SMsort[i],i] <- tanvec[SMsort[i],3]
        }
    }
    U <- weights*U
    surfOnly <- surface
    freeOnly <- free
    if (!is.null(surface))
        surfOnly <- unique(sort(surface))
    if (!is.null(free))
        freeOnly <- unique(sort(free))
    ## remove fix columns
    
    Ured0 <- as(U[,1:m],"sparseMatrix")
    
    if (!is.null(surface) || !is.null(free)) {
        Ured1 <- as(U[,(m+1):(2*m)],"sparseMatrix")
        smsurffree <- which(! SMvector %in% c(surfOnly,freeOnly))
        if (length(smsurffree))
            Ured1 <- Ured1[,-smsurffree]
        Ured0 <- cBind(Ured0,Ured1)
        
    }
    if (!is.null(free)) {
        Ured1 <- as(U[,(2*m+1):(3*m)],"sparseMatrix")
        smsurffree <- which(! SMvector %in% c(freeOnly))
        if (length(smsurffree))
            Ured1 <- Ured1[,-smsurffree]
        Ured0 <- cBind(Ured0,Ured1)
        
    }
    
    return(list(tanvec=tanvec,SMvector=SMvector,Gamma0=Gamma0,U=Ured0))             
}     

