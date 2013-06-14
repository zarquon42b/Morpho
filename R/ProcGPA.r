ProcGPA<-function(dat.array,tol=1e-5,scale=TRUE,CSinit=FALSE,silent=FALSE,weights=NULL,centerweight=FALSE)
{
  if (!is.null(weights))
    weights <- weights/sum(weights)

  t0<-Sys.time()
  x<-dat.array
  p1<-1e10
  p2<-p1	
  n<-dim(dat.array)[3]
  k<-dim(dat.array)[1]
  m<-dim(dat.array)[2]
  x1<-gdif(dat.array)
                                       
  arr.list<-list(0)	
###rotation step ####
  for ( i in 1:n)
    {
      arr.list[[i]]<-list(x[,,i],1)
    }
  
  if (CSinit)
    {
      arr.list<-lapply(arr.list,function(x){x[[1]]<-apply(x[[1]],2,scale,scale=F);x[[1]]<-x[[1]]/sqrt(sum(x[[1]]^2));return(list(x[[1]],x[[2]]))})
    }
  
  else 
    {
      arr.list<-lapply(arr.list,function(x){x[[1]]<-apply(x[[1]],2,scale,scale=F);return(list(x[[1]],x[[2]]))})
    }
  
  mshape<-x[,,1]
   if (centerweight && !is.null(weights))
          {
            mcent <- apply(mshape*weights,2,sum)           
            mshape<-scale(mshape,scale=F,center=mcent)
          }
### align mean by principal axes ###	
  rotms<-eigen(crossprod(mshape))$vectors
  if (det(rotms) < 0)
    {
      rotms[,1]<-rotms[,1]*-1
    }
  mshape<-mshape%*%rotms
  
  while (p1 > tol)
    {
      mshape_old<-mshape
      
### rotation of all configs on current consensus ###		
      arr.list<-lapply(arr.list,function(x){x[[1]]<-rot.proc(x[[1]],x=mshape,scale=F,weights=weights,centerweight=centerweight);return(list(x[[1]],x[[2]]))})
      
      for( i in 1:n)
        {
          x[,,i]<-arr.list[[i]][[1]]
        }
      
      x2<-gdif(x)
      p1<-abs(x1-x2)
      x1<-x2
    }
  
### scale/rotate step ###	
  if (scale)
    {      
      for ( i in 1:n)
        {
          arr.list[[i]]<-list(x[,,i],1)
        }
	
      while (p2 > tol)
        {
          for( i in 1:n)
            { if (!is.null(weights))
                x[,,i]<-arr.list[[i]][[1]]*weights
            else
              x[,,i]<-arr.list[[i]][[1]]
            }
          eigc<-scaleproc(x)
          
          for ( i in 1:n)	
            {
              arr.list[[i]][[2]]<-eigc[i]
            }
          
          arr.list<-lapply(arr.list,function(x){x[[1]]<-x[[1]]*x[[2]];return(list(x[[1]],x[[2]]))})         
                    
### rotation of all configs on current consensus ###		
          arr.list<-lapply(arr.list,function(x){x[[1]]<-rot.proc(x[[1]],x=mshape,scale=F,weights=weights,centerweight=centerweight);return(list(x[[1]],x[[2]]))})
          		
### scale step ####
          
          for( i in 1:n)
            {
              x[,,i]<-arr.list[[i]][[1]]
            }	
          x2<-gdif(x)
          p2<-abs(x1-x2)
          x1<-x2
        }
    }
  mshape<-apply(x,c(1,2),mean)
  if (CSinit)
    {
      msize<-cSize(mshape)
      mshape<-mshape/msize
      if (scale)
        {
          x<-x/msize
        }
    }
  t1<-Sys.time()
  if (!silent)
    {
      cat(paste("in... ",format(t1-t0)[[1]],"\n"))
    }	
  return(list(rotated=x,mshape=mshape))
  
}	
