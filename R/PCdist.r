#' correlation between a reduced space and the original space
#' 
#' Calculates the correlation between distances in a reduced space and the
#' original space
#' 
#' 
#' @param PCs m x k matrix of Principal Components where m is the k is the
#' number of PCs.
#' @param PCscores n x m matrix of Principal Component scores where n is the
#' number of observations.
#' @param x integer: increment for every x-th PC the subspace to fullspace
#' correlation will be calculated.
#' @param plot.type "b"=barplot of correlation values, "s"=line between
#' correlation values.
#' @return a vector of R-squared values between subspace and fullspace
#' distances and a barplot depicting the correlations belonging to the
#' subspace.
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' library(shapes)
#' a <- procSym(gorf.dat)
#' PCdist(a$PCs, a$PCscores, x = 2)
#' 
#' 
#' @export
PCdist <- function(PCs,PCscores,x=5,plot.type="b")
{

    k <- dim(PCs)[2]
    rest <- k%%x
    mod <- floor(k/x)
    if (rest==0)
        {bar <- mod}
    else
        {bar <- mod+1}
    alld <- t(PCs%*%t(PCscores))
    alldist <- dist(alld)
    dist.list <- list(numeric(0))
    cor.vec <- c(numeric(0))

    if (rest==0)
      { nam <- c(seq(from=x,by=x,length.out=mod))
        for (i in 1:bar)
        {
          part <- t(PCs[,1:(i*x)]%*%t(PCscores[,1:(i*x)]))
          dist.list[[i]] <- dist(part)
          cor.vec[i] <- (cor(dist.list[[i]],alldist))^2
        }
      }
    else
      { nam <- c(seq(from=x,by=x,length.out=mod),k)
        for (i in 1:(bar-1))
        {
          part <- t(PCs[,1:(i*x)]%*%t(PCscores[,1:(i*x)]))
          dist.list[[i]] <- dist(part)
          cor.vec[i] <- (cor(dist.list[[i]],alldist))^2
        }
        dist.list[[bar]] <- alldist
        cor.vec[bar] <- 1
      }
      if (plot.type=="b")
      
      {barplot(cor.vec,xlab="PCs",ylab="correlation between reduced space and full space (R-squared)",names.arg=nam,xpd=FALSE)}
      
      else 
            
      { plot(nam,cor.vec,xlab="PCs",ylab="R-squared",main="correlation between reduced space and full space",cex=0.7)
        text(nam,cor.vec,labels=nam,cex=0.7, font=4,pos=3)
        lines(nam,cor.vec)}
      
      
      return(cor.vec)
}

    

      

    
