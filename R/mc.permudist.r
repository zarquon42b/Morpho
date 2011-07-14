mc.permudist <- function(data,groups,rounds,which=1:2)
  {
    require(foreach)
    require(doMC)
    registerDoMC()
   
    N <- data
    if (is.vector(N))
      {N <- as.matrix(N)
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
    #print(groups)
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
#print(l2)
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
       # print(shake)
        group.tmp[[1]] <- shake[1:l1]
        group.tmp[[2]] <- shake[(l1+1):lshake]
      # print(group.tmp)
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

    return(list(permudist=dists,dist=dist))
  }
