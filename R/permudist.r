mc.permudist <- permudist <- function(data,groups,rounds=1000,which=1:2,mc.cores = detectCores())
  {
     if(.Platform$OS.type == "windows")
      registerDoParallel(makeCluster(1),cores=1)
    else
      registerDoParallel(cores=mc.cores)### register parallel backend
   
### configure grouping ####
    N <- data
    if (is.vector(N))
      {N <- as.matrix(N)
     }
    if (dim(N)[2] == 3)
      {N <- vecx(N)
     }
    if (!is.factor(groups))
      {groups<-as.factor(groups)
     }
    if (is.factor(groups))
      {
        lev<-levels(groups)
        levn<-length(lev)
        group<-list()
        count<-1
        groupcheck<-0
        for (i in 1:levn)
          {	tmp0<-which(groups==lev[i])	
                if (length(tmp0) != 0)
                  {			
                    group[[count]]<-tmp0
                    count<-count+1
                    groupcheck[count]<-i
                  }
              }
        groups<-group
      }
### end configure grouping ####

    b <- groups
    n <- dim(N)[1]
    
    l <- dim(N)[2]
    l1 <- length(groups[[which[1]]])
    l2 <- length(groups[[which[2]]])
    if (dim(N)[2] == 1)
      {
        mean1 <- mean(N[groups[[which[1]]]])
        mean2 <- mean(N[groups[[which[2]]]])
      }
    else
      {
        mean1 <- apply(N[groups[[which[1]]],],2,mean)
        mean2 <- apply(N[groups[[which[2]]],],2,mean)
      }
    dist <- sqrt(sum((mean1-mean2)^2))
    dists <- NULL
    
    permu <- function(x)
      {
        group.tmp <- list()
        shake <- sample(unlist(groups[which]))
        lshake <- length(shake)
        
        group.tmp[[1]] <- shake[1:l1]
        group.tmp[[2]] <- shake[(l1+1):lshake]
     
        if (dim(N)[2] == 1)
          {
            mean1.tmp <- mean(N[group.tmp[[1]]])
            mean2.tmp <- mean(N[group.tmp[[2]]])
          }
        else
          {
            mean1.tmp <- apply(N[group.tmp[[1]],],2,mean)
            mean2.tmp <- apply(N[group.tmp[[2]],],2,mean)
          }
        disto <- sqrt(sum((mean1.tmp-mean2.tmp)^2))
        return(disto)
      }
    dists <- foreach(i= 1:rounds,.combine=c) %dopar%
    permu(i)
    p.value <- length(which(dists >= dist))
    if (p.value > 0)
      {
        p.value <- p.value/rounds
        names(p.value) <- "p-value"
      }
    else
      {
        p.value <- 1/rounds
         names(p.value) <- "p-value <"
      }
    
    return(list(permudist=dists,dist=dist,p.value=p.value))
  }
