rot.proc<-function(x,y,scale=TRUE)
	{
	XY<-crossprod(x,y)	
	sv1<-svd(XY)
  	#dd<-diag(sign(sv1$d))
	gamm<-tcrossprod(sv1$v,sv1$u)
	#gamm<-(sv1$v)%*%gamm
	del<-sv1$d
  
  if (scale == TRUE)
    {	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
  	bet<-sum(del)/ctrace(y)
	yrot<-bet*y%*%gamm}
  else
    	{yrot<-y%*%gamm}
	return(yrot)
	}
