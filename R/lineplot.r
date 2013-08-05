lineplot  <- function(x,point,col=1,lwd=1,line_antialias = FALSE,add=TRUE)
{
  
  if (dim(x)[2] == 3)
    {
      if (is.list(point)==TRUE)
        {
          for (i in 1:length(point))
            {
              lines3d(x[point[[i]],1],x[point[[i]],2],x[point[[i]],3],col=col,lwd=lwd,line_antialias = line_antialias)
           }
        }
      else                                                
        {
          lines3d(x[point,1],x[point,2],x[point,3],col=col,lwd=lwd,line_antialias = line_antialias)
        }
    }
  else
    {
      if (!add)
        {
          plot(x,asp=1,cex=0,xlab="x-coordinate",ylab="y-coordinate")
        }
      if (is.list(point)==TRUE)
        {
          for (i in 1:length(point))
            {
              lines(x[point[[i]],1],x[point[[i]],2],col=col,lwd=lwd)
            }
        }
      else                                                
        {
          lines(x[point,1],x[point,2],col=col,lwd=lwd)
      }
    }
}
