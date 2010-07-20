lineplot3<-function(x,point,col=1,lwd=1,line_antialias = FALSE)
{     

  if (is.list(point)==TRUE)
    {for (i in 1:length(point))
      {
      lines3d(x[point[[i]],1],x[point[[i]],2],x[point[[i]],3],col=col,lwd=lwd,line_antialias = line_antialias)
      }}
     else                                                
     {lines3d(x[point,1],x[point,2],x[point,3],col=col,lwd=lwd,line_antialias = line_antialias)
     }
}                    