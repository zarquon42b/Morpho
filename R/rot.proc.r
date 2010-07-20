rot.proc<-function(x,y,scale=TRUE)
	{
	XY<-crossprod(x,y)	
	sv1<-svd(XY)
  	gamm<-(sv1$v)%*%t(sv1$u)
	del<-sv1$d
  
  if (scale == TRUE)
    {	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
  	bet<-sum(del)/ctrace(y)
	bet<-sum(del)/ctrace(y)
	yrot<-bet*y%*%gamm}
  else
    	{yrot<-y%*%gamm}
	return(yrot)
	}
