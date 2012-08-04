rot.proc<-function(x,y,scale=TRUE,weights=NULL,centerweight=FALSE)
  {
   
     if (centerweight && !is.null(weights))
          {
            xcent <- apply(x*weights,2,sum)
            ycent <- apply(y*weights,2,sum)
            x<-scale(x,scale=F,center=xcent)
            y<-scale(y,scale=F,center=ycent)
          }
### rotates 2 already centred matrices onto each other#
     if (!is.null(weights))
          {
            Dn <- diag(weights)
            X1 <- Dn%*%x
            Y1 <- Dn%*%y
            XY <- crossprod(X1,Y1)
          }
     else
       XY<-crossprod(x,y)	
     sv1<-svd(XY)
                                        #dd<-diag(sign(sv1$d))
    gamm<-tcrossprod(sv1$v,sv1$u)
                                        #gamm<-(sv1$v)%*%gamm
    del<-sv1$d
    
    if (scale == TRUE)
      {	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
         if (!is.null(weights))
           bet <- sum(del)/ctrace(Y1)
         else
           bet<-sum(del)/ctrace(y)
	yrot<-bet*y%*%gamm}
    else
      {yrot<-y%*%gamm}
    return(yrot)
  }
