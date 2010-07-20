rotonto<-function(x,y,scaling=FALSE,signref=TRUE)
{ 
  	m<-dim(x)[2]
  	X<-apply(x,2,scale,scale=F)
  	Y<-apply(y,2,scale,scale=F)
  	XY<-crossprod(X,Y)
  	sv1<-svd(XY)
  	gamm<-(sv1$v)%*%t(sv1$u)
  	if(sign(det(gamm))<1 && signref ==T)
  		{cat("reflection involved\n")
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
  	yrot<-t(apply(yrot,1,function(x){x+trans}))
  
  	return(list(yrot=yrot,Y=Y,X=X,trans=trans,transy=transy,gamm=gamm,bet=bet))
}
