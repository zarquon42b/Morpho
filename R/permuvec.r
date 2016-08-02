#' perfom permutation testing on angles and distances between subgroups of two
#' major groups.
#' 
#' 
#' perform permutation test on length and angle of the vectors connecting the
#' subgroup means of two groups: e.g. compare if length and angle between sex
#' related differences in two populations differ significantly.
#' 
#' This function calculates means of all four subgroups and compares the
#' residual vectors of the major grouping variables by angle and distance.
#' 
#' @param data array or matrix containing data.
#' @param groups factors of firs two grouping variables.
#' @param subgroups factors of the subgrouping.
#' @param rounds number of requested permutation rounds
#' @param scale if TRUE: data will be scaled by pooled within group covarivance
#' matrix. Otherwise Euclidean distance will be used for calculating distances.
#' @param tol threshold for inverting covariance matrix.
#' @param mc.cores integer: determines how many cores to use for the
#' computation. The default is autodetect. But in case, it doesn't work as
#' expected cores can be set manually.Parallel processing is disabled on
#' Windows due to occasional errors.
#' @return
#' \item{angle }{angle between the vectors of the subgroups means}
#' \item{dist }{distances between subgroups}
#' \item{meanvec }{matrix containing the means of all four subgroups}
#' \item{permutangles }{vector containing angles (in radians) from random permutation}
#' \item{permudists }{vector containing distances from random permutation}
#' \item{p.angle }{p-value of angle between residual vectors}
#' \item{p.dist }{p-value of length difference between residual vectors}
#' \item{subdist }{length of residual vectors connecting the subgroups}
#' means.
#' 
#' @examples
#' 
#' data(boneData)
#' proc <- procSym(boneLM)
#' pop <- name2factor(boneLM,which=3)
#' sex <- name2factor(boneLM,which=4)
#' ## use non scaled distances by setting \code{scale = FALSE}
#' ## and only use first 10 PCs
#' perm <- permuvec(proc$PCscores[,1:10], groups=pop, subgroups=sex,
#'                  scale=FALSE, rounds=100, mc.cores=2)
#' 
#' 
#' ## visualize if the amount of sexual dimorphism differs between
#' # (lenghts of vectors connecting population specific sex's averages)
#' # differs between European and Chines
#' hist(perm$permudist, xlim=c(0,0.1),main="measured vs. random distances",
#'      xlab="distances")
#' points(perm$dist,10,col=2,pch=19)#actual distance
#' text(perm$dist,15,label=paste("actual distance\n
#'      (p=",perm$p.dist,")"))
#' ## not significant!!
#' 
#' ## visualize if the direction of sexual dimorphism
#' # (angle between vectors connecting population specific sex's averages)
#' # differs between European and Chines
#' hist(perm$permutangles, main="measured vs. random angles",
#'      xlab="angles")
#' points(perm$angle,10,col=2,pch=19)#actual distance
#' text(perm$angle,15,label=paste("actual distance\n
#'     (p=",perm$p.angle,")"))
#' ## also non-significant
#' 
#' @export
permuvec <- function(data,groups,subgroups=NULL,rounds=10000,scale=TRUE,tol=1e-10,mc.cores=parallel::detectCores())
{
  win <- FALSE
  if(.Platform$OS.type == "windows")
    win <- TRUE
  else
    registerDoParallel(cores=mc.cores)### register parallel backend
  
### define groups ####
  rawgroup <- groups	
  lev <- NULL	
  if (is.character(groups) || is.numeric(groups))
      groups <- as.factor(groups)
  if (is.factor(groups)) {
      lev <- levels(groups)
      levn <- length(lev)
      group <- list()
      count <- 1
      groupcheck <- 0
      for (i in 1:levn) {
          tmp0 <- which(groups==lev[i])	
          if (length(tmp0) != 0) {			
              group[[count]] <- tmp0
              count <- count+1
              groupcheck[count] <- i
          }
      }
      groups <- group
  }
  levsub <- NULL	
  if (is.character(subgroups) || is.numeric(subgroups))
      subgroups <- factor(subgroups)
  
  if (is.factor(subgroups)) {
      levsub <- levels(subgroups)
      levnsub <- length(levsub)
      subgroup <- list()
      count <- 1
      for (i in 1:levnsub) {
          tmp0 <- which(subgroups==levsub[i])
          if (length(tmp0) != 0) {
              subgroup[[i]] <- tmp0
              count <- count+1
          }
      }
      subgroups <- subgroup
  }
  b <- groups
  N <- data
  ng <- length(groups)
  nsub <- length(subgroups)
  meanlist <- list()
  
### prepare data if data is an array ###
  
  if (length(dim(N)) == 3) {
      n <- dim(N)[3]
      k <- dim(N)[1]
      m <- dim(N)[2]
      l <- k * m
      
      if (length(unlist(groups)) != n)
          warning("group affinity and sample size not corresponding!")
      
      nwg <- c(rep(0, ng))
      for (i in 1:ng) 
          {nwg[i] <- length(b[[i]])
       }
      
      B <- matrix(0, n, m * k)
      for (i in 1:n) 
          B[i, ] <- as.vector(N[, , i])
  } else {
      n <- dim(N)[1]
      l <- dim(N)[2]
      if (length(unlist(groups)) != n)
          warning("group affinity and sample size not corresponding!")
      
      ng <- length(groups)
      nwg <- c(rep(0, ng))
      for (i in 1:ng)
          nwg[i] <- length(b[[i]])
      
      B <- as.matrix(N)
  }
  nws <- c(rep(0, nsub))
  for (i in 1:nsub)
      nws[i] <- length(subgroups[[i]])
  Gmeans <- matrix(0, ng, l) ### calculate mean of subgroup means for all groups ###
  for (i in 1:ng) {
      for (j in 1:nsub) {
          tmp <- subgroups[[j]][which(subgroups[[j]] %in% groups[[i]])]
          Gmeans[i,] <- Gmeans[i,]+colMeans(B[tmp,,drop=FALSE])
      }
      Gmeans[i,] <- Gmeans[i,]/nsub
  }
  
### create empty subgroup mean matrices 
  meanvec <- matrix(NA,ng,l)
  if (is.factor (rawgroup ) || is.character(rawgroup))
      rownames(meanvec) <- levels(rawgroup)[groupcheck]
  for (i in 1:ng)
      meanlist[[i]] <- matrix(NA,nsub,l)
   
### correct for groupmeans ###
  for (i in 1:ng) {
      gn <- length(groups[[i]])
      delt <- matrix(Gmeans[i,],gn,l,byrow=T)
      B[groups[[i]],] <- B[groups[[i]],]-delt
  }
  
### calc subgroup means, residual vectors and pooled within group variance ###
  covW <- 0	
  for (i in 1:ng) {	
      for (j in 1:nsub) {
          tmp <- subgroups[[j]][which(subgroups[[j]] %in% groups[[i]])]
          meanlist[[i]][j,] <- colMeans(B[tmp,,drop=FALSE])
### calc within subgroups Sum of Squares
          if (scale)
              covW <- covW+cov(scale(B[tmp,,drop=FALSE], scale=F))*(length(tmp)-1)
      }
      
### calc pooled groupspecific within subgroups covariance matrix and overall variance ###
      meanvec[i,] <- (meanlist[[i]][1,]-meanlist[[i]][2,])      
  }
  covW <- covW/(n-(ng*nsub))
  if (!scale)
      covW <- diag(rep(1,dim(B)[2]))
  mahadist <- NULL
### invert covariance matrix
  coinv <- armaGinv(covW,tol=tol)
                                        # print(dim(coinv))
  for (i in 1:ng) ## calc Mahalanobisdistance ### 
    mahadist[i] <- sqrt(meanvec[i,]%*%coinv%*%meanvec[i,])
      
### calc angle compare vector lengths ###
  
  disto <- abs(mahadist[1]-mahadist[2])
  out <- angle.calc(meanvec[1,],meanvec[2,])
  
### permutate over groups ###	
  
  permuta <- function(x)
    {
      tmplist <- list()
      for (i in 1:ng)
          tmplist[[i]] <- matrix(NA,nsub,l)
      meanvectmp <- matrix(NA,ng,l)
      Btmp <- B
      b1 <- list(numeric(0))
      shake <- sample(1:n)
      #Gmeans1 <- matrix(0, ng, l)
      l1 <- 0
      for (j in 1:ng) {
          b1[[j]] <- c(shake[(l1 + 1):(l1 + (length(b[[j]])))])
          l1 <- l1 + length(b[[j]])
      }
      for (i in 1:ng) {
          for (j in 1:nsub) {
              tmp <- subgroups[[j]][which(subgroups[[j]] %in% b1[[i]])]
              tmplist[[i]][j,] <- colMeans(Btmp[tmp,,drop=FALSE])
          }
          meanvectmp[i,] <- (tmplist[[i]][1,]-tmplist[[i]][2,])
      }
      
      mahadist0 <- NULL
      for (i in 1:ng) ## calc Mahalanobisdistance ### 
          mahadist0[i] <- sqrt(meanvectmp[i,]%*%coinv%*%meanvectmp[i,])        
        
      dist <- abs(mahadist0[1]-mahadist0[2])
      
      return(c(angle.calc(meanvectmp[1,],meanvectmp[2,]),dist))
  }
  if (win)
      tt <- foreach(i= 1:rounds,.combine=c) %do% permuta(i)
  else
      tt <- foreach(i= 1:rounds,.combine=c) %dopar% permuta(i)
                                        # mclapply(alist,permuta)
  uns <- unlist(tt)
  angs <- (1:rounds)*2-1
  dists <- uns[2*(1:rounds)]
  sortdist <- sort(dists)
  
### calc probabilities ####
  if (max(sortdist) < disto) {
      probadist <- 1/rounds
    } else {
        marg <- min(which(sortdist >= disto))
        probadist <- (rounds-marg+1)/rounds
    }      
  sortang <- sort(uns[angs])
  if (max(sortang) < out) {
      proba <- 1/rounds
  } else {
      marg <- min(which(sortang >= out))
      proba <- (rounds-marg+1)/rounds
  }      
  
  return(list(angle=out,dist=disto,meanvec=meanvec,permutangles=sortang,permudists=sortdist,p.angle=proba,p.dist=probadist,subdist=mahadist))
}

