difplot.proc<-function(A,color=4,lwd=1,away.fac=1.1,lcol=2,sz=0,output=TRUE)
{   rgl.clear()
    rgl.bg(color = "white")
    n<-dim(A$rotated)[3]
    disti<-data.frame(c(1:n),A$rho)
    disti.sort<-disti[order(disti[,2],decreasing=T),]
    colnames(disti.sort)[1]<-"# in array"
    rownames(disti.sort)<-c(1:n)
    outlier<-0
    t1 <- 1
    lo<-1
    answer1 <- "n"
    while (t1 <= n )
    {
      if (answer1 %in% c("n","N","y","Y") == T )
      {
      difplot(A$mshape,A$rotated[,,disti.sort[t1,1]],ID=disti.sort[t1,1],color=color,lwd=1,lcol=lcol,sz=0,rgl.new=FALSE,away.fac=away.fac)
      
      cat(paste("outlier #",t1,": ",disti.sort[t1,1],"     procrustes dist. to mean: ",disti.sort[t1,2],"\n",sep=""))
      
      answer <- substr(readline(" add to outlierlist (y/N)?  "), 1L,1L)
       if (answer == "y" || answer == "Y")
         {outlier[lo]<-disti.sort[t1,1]
          lo<-lo+1
          }
      # else if (answer == "n" || answer == "N")
          
      }
      
      answer1 <- substr(readline(" wanna go on? (y/N)?  "), 1L,1L)
      
      if (answer1 == "y" || answer1 == "Y")
         {
           rgl.clear()
           rgl.bg(color = "white")
            t1<-t1+1          
         }
      else if (answer1 == "n" || answer1 == "N")
          {break}
      
      else {cat("thy answer be yes or no!! \n")
            answer1 <- "q"
            }
     
    }
    
    
    
    
    if (output == TRUE)
    {    return(list(outlier=outlier,dist.sort=disti.sort))}
    
}
