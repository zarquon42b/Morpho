rotonto<-function(x,y,scaling=FALSE,signref=TRUE,reflection=TRUE)
{ 	reflect=0
  	m<-dim(x)[2]
  	X<-apply(x,2,scale,scale=F)
  	Y<-apply(y,2,scale,scale=F)
  	XY<-crossprod(X,Y)
  	sv1<-svd(XY)
        
	#dd<-diag(sign(sv1$d))
	gamm<-tcrossprod(sv1$v,sv1$u)
	#gamm<-(sv1$v)%*%gamm
  	
  	if(sign(det(gamm))<1)
	  {	reflect<-1
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
                
                  
	  
  	trans<-x[1,]-X[1,]
  	transy<-y[1,]-Y[1,]
  	#yrot<-Y%*%gamm
  	sig<-sign(det(XY))
  	del<-sv1$d
   #del[m]<-sig*abs(del[m])
  	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
  	bet<-sum(del)/ctrace(Y)
  
  	if (scaling == TRUE)
    		{yrot<-bet*Y%*%gamm
		}
  	else
    		{yrot<-Y%*%gamm
		}
	Y<-yrot  	
	yrot<-t(apply(yrot,1,function(x){x+trans}))
  
  	return(list(yrot=yrot,Y=Y,X=X,trans=trans,transy=transy,gamm=gamm,bet=bet,reflect=reflect))
}
