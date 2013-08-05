mcNNindex <- function(target,query,cores=detectCores(),k=k,...)
    {
        if(.Platform$OS.type == "windows")
            cores <- 1
        if (nrow (query) < 1000)
            cores <- 1
        out <- NULL
        mclist <- list()
        nx <- dim(query)[1]
        iter <- floor(nx/cores)
        if (cores > 1)
            {
                for (i in 1:(cores-1))
                    mclist[[i]] <- query[(1:iter)+((i-1)*iter),]
                
                mclist[[cores]] <- query[-c(1:((cores-1)*iter)),]
            }
        else
            mclist[[1]] <- query
        tmpfun <- function(x,...)
            {
                ##tmp0 <- nn2(target,x,k=k,searchtype="priority",...)$nn.idx
                tmp0 <- ann(ref=target, target=x, k=k, search.type="priority",verbose=FALSE)$knnIndexDist[,1:k] ## ann function from package yaImpute
                return(tmp0)
            }
        tmp <- mclapply(mclist,tmpfun,mc.cores=cores)
        for (i in 1:cores)
            out <- rbind(out,tmp[[i]])
            
        return(out)

    }
