.calcTang_U_s <- function(datamatrix,normalmatrix=NULL,SMvector,outlines=NULL,surface=NULL,free=NULL,deselect=FALSE)
{
    
    dims <- dim(datamatrix)[2]
    k <- dim(datamatrix)[1]
        
    if (deselect==TRUE)
        SMvector <- c(1:k)[-SMvector]
    
    SMvector <- unique(SMvector)
    m <- length(SMvector)
    if ( !is.null(free)) {
        udims <- c(dims*k,m*3)
        tanvec <- matrix(0,k,dims*3)
        #U <- matrix(0,dims*k,m*3)
        type <- 2
    } else if(!is.null(surface)) {
        udims <- c(dims*k,m*2)
        tanvec <- matrix(0,k,dims*2)
        #U <- matrix(0,dims*k,m*2)
        type <- 1
    } else {
        udims <- c(dims*k,m)
        tanvec <- matrix(0,k,dims)
        #U <- matrix(0,dims*k,m)
        type <- 0
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
            tanp <- tangentPlane(normalmatrix[temp[i],])
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
    xinfo <- .Call("tweakU",tanvec,m, type,SMsort)
    U <- sparseMatrix(i=xinfo$rows,j=xinfo$cols+1, x=xinfo$x,dims=udims)
    #U <- weights*U
    
    
    outOnly <- outlines
    surfOnly <- surface
    freeOnly <- free
    if (!is.null(surface))
        surfOnly <- unique(sort(surface))
    if (!is.null(free))
        freeOnly <- unique(sort(free))
    if(!is.null(outlines))
        outOnly <- unique(sort(unlist(outlines)))
    
    allsurf <- c(outOnly,surfOnly,freeOnly)
    
    if (length(which(!SMvector %in% allsurf)))
        stop("all semi-landmarks must to be tagged as outlines or surfaces")
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
    
    return(list(SMvector=SMvector,Gamma0=Gamma0,U=Ured0))             
}     

