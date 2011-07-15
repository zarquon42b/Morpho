plot.curvature <- function(curve,text=FALSE,color="rainbow")
  {
    limits <- range(curve$id)
    if (color !="rainbow")
      {points3d(curve$vertices,col=color)
       if (text || color == "rainbow")
         {
           for (i in curve$id)
             {
               if (color == "rainbow")
                 {
                   color=i+1
                   points3d(curve$vertices[which(curve$id==i),],col=color)
                 }
               
               
               if (text)
                 {
                   text3d(apply(curve$vertices[which(curve$id==i),],2,mean),texts=paste(i),col=color)
                 }
             }
         }
     }
  }
