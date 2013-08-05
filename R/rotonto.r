rotonto <- function(x,y,scale=FALSE,signref=TRUE,reflection=TRUE,weights=NULL,centerweight=FALSE)
{ 	reflect=0
  	m <- dim(x)[2]
        if (!is.null(weights))
          weights <- weights/sum(weights)
       
        X <- apply(x,2,scale,scale=F)
  	Y <- apply(y,2,scale,scale=F)
         if (centerweight && !is.null(weights))
          {
            
            xcent <- apply(X*weights,2,sum)
            ycent <- apply(Y*weights,2,sum)
            X <- scale(X,scale=F,center=xcent)
            Y <- scale(Y,scale=F,center=ycent)
          }
        if (!is.null(weights))
          {
            Dn <- diag(weights)
            X1 <- Dn%*%X
            Y1 <- Dn%*%Y
            XY <- crossprod(X1,Y1)
          }
        else
          XY <- crossprod(X,Y)

        sv1 <- svd(XY)
        
	#dd <- diag(sign(sv1$d))
	gamm <- tcrossprod(sv1$v,sv1$u)
	#gamm <- (sv1$v)%*%gamm
  	
  	if(sign(det(gamm))<1)
	  {	reflect <- 1
		if (signref && reflection)
                  {cat("reflection involved\n")
                 }
                if (!reflection)
                  {
                    u <- sv1$u
                    v <- sv1$v
                    chk1 <- Re(prod(eigen(v)$values))
                    chk2 <- Re(prod(eigen(u)$values))
                    if ((chk1 < 0) && (chk2 > 0))
                      {
                        v[, dim(v)[2]] <- v[, dim(v)[2]] * (-1)
                        gamm <- v %*% t(u)
                      }
                    if ((chk2 < 0) && (chk1 > 0))
                      {
                        u[, dim(u)[2]] <- u[, dim(u)[2]] * (-1)
                        gamm <- v %*% t(u)
                      }
                  }
              }
                
                  
	  
  	trans <- x[1,]-X[1,]
  	transy <- y[1,]-Y[1,]
  	#yrot <- Y%*%gamm
  	sig <- sign(det(XY))
  	del <- sv1$d
   #del[m] <- sig*abs(del[m])
  	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
  
  
  	if (scale == TRUE)
    		{
                  if (!is.null(weights))
                    bet <- sum(del)/ctrace(Y1)
                  else
                    bet <- sum(del)/ctrace(Y)
                  yrot <- bet*Y%*%gamm
		}
  	else
    		{
                  bet <- 1
                  yrot <- Y%*%gamm
		}
	Y <- yrot  	
	yrot <- t(apply(yrot,1,function(x){x+trans}))
  
  	return(list(yrot=yrot,Y=Y,X=X,trans=trans,transy=transy,gamm=gamm,bet=bet,reflect=reflect))
}
rotreverse <- function(mat,rot)
  {
    transfun <- function(x,trans)
      {
        x <- x+trans
      }
    
    out <- t(apply(mat,1,transfun,-rot$trans))
    out <- t(apply((out%*%t(rot$gamm))*1/rot$bet,1,transfun,rot$transy))
    return(out)
  }
