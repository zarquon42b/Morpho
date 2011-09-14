mc.pls2B <- function(y,x,tol=1e-12,rounds=0)
  {
    require(foreach)
    require(doMC)
    registerDoMC()
    
    
    if (length(dim(x)) == 3)
      {
        x <- vecx(x)
      }
    
    xdim <- dim(x)
    ydim <- dim(y)
   
      
    
    cova <- cov(cbind(x,y))
    svd.cova <- svd(cova[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])

    svs <- svd.cova$d
    svs <- svs/sum(svs)
    svs <- svs[which(svs > 0.001)]

    covas <- svs*100
    l.covas <- length(covas)
    z1 <- x%*%svd.cova$u[,1:l.covas] #pls scores of x
    z2 <-  y%*%svd.cova$v[,1:l.covas] #pls scores of y
    
### calculate correlations between pls scores
    cors <- 0
    for(i in 1:length(covas))
      {cors[i] <- cor(z1[,i],z2[,i])
     }



### Permutation testing
   
      
    permupls <- function(i)
      {
        x.sample <- sample(1:xdim[1])
        y.sample <- sample(x.sample)
       
        cova.tmp <- cov(cbind(x[x.sample,],y[y.sample,]))
        svd.cova.tmp <- svd(cova.tmp[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])
        svs.tmp <- svd.cova.tmp$d
        return(svs.tmp[1:l.covas])
      }
    p.values <- rep(NA,l.covas)
    if (rounds > 0)
      {
        permuscores <- foreach(i = 1:rounds, .combine = cbind) %dopar% permupls(i)
        
        p.val <- function(x,rand.x)
          {
            p.value <- length(which(rand.x >= x))
                       
            if (p.value > 0)
              {
                p.value <- p.value/rounds
              }
            else
              {p.value <- 1/rounds}
            gc()
            return(p.value)
          }
            
            
            for (i in 1:l.covas)
              {
                p.values[i] <- p.val(svd.cova$d[i],permuscores[i,])
                
               
               }
          }
        ### create covariance table
    
    Cova <- data.frame(svd.cova$d[1:l.covas],covas,cors,p.values)
    colnames(Cova) <- c("singular value","% total covar.","Corr. coefficient", "p-value")
    out <- list(svd=svd.cova,z2=z2,z1=z1,cor=cors,CoVar=Cova)
    return(out)
  }
