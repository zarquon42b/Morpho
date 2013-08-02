fixLMtps<-function(data,comp=3,weight=TRUE)
{
  n<-dim(data)[3]
  k<-dim(data)[1]
  m<-dim(data)[2]
  checklist<-list()
  checkvec<-rep(0,n)
  out<-data
  ## check for missing landmarks ###
  for (i in 1:n)
    {count<-0
     found<-FALSE
     checklist[[i]]<-NA
     for (j in 1:k)
       {
         if (NA%in%data[j,,i])
           {count<-count+1
            checklist[[i]][count]<-j
            checkvec[i]<-1
          }
       }
   }
  ## calc mean of complete configs ###
  check<-which(checkvec==1)
  data.c<-data[,,-check]
  ## check if there are enough configs to use weighting option
  if (length(dim(data.c)) < 3)
    {
      if (is.matrix(data.c))
        {
          data.c <- array(data.c,dim=c(dim(data.c),1))
          ngood <- 1
        }
      else
        stop("there is no complete configuration to use")
    }
  else
    ngood <- dim(data.c)[3]
  if (ngood < comp)
    {
      if (ngood == 0)
        stop("no complete configuration found")
      else
        { comp <- ngood
          if (weight)
            warning(paste("only",ngood,"configurations found. comp is set to",ngood,"\n"))
        }
    }
  
  ## rotate incomplete data onto mean ###
  lmsdat<-data
  if (ngood > 1)
    {
      proc.c<-ProcGPA(data.c,silent = TRUE)
      mean0<-proc.c$mshape
    }
  else
    mean0 <- data.c[,,1]
  for (i in 1:length(check))
    {
      miss<-checklist[[check[i]]]

      if (weight && ngood > 1) ### calculate weighted estimates of missing data ###
        {
          ## rotate incomplete data onto mean ###
          rotmiss <- rotonto(mean0[-miss,],data[-miss,,check[i]],scale=TRUE)$yrot
          allrot <- bindArr(rotmiss,proc.c$rotated[-miss,,], along=3)
          ## calculate weights according to procrustes distance ###			
          wcalc <- proc.weight(allrot,comp,1,report=FALSE)
          lms <- proc.c$rotated[,,wcalc$data$nr-1]
          lm.est <- matrix(0,dim(data)[1],m)
          
          for (j in 1:comp)
            {lm.est<-lm.est+lms[,,j]*wcalc$data$weight[j]
           }
          tpsout<-tps3d(lm.est,lm.est[-miss,],data[-miss,,check[i]])
        }
      else
        {
          tpsout<-tps3d(mean0,mean0[-miss,],data[-miss,,check[i]])
        }
      out[,,check[i]]<-tpsout
    }
  return(list(out=out,mshape=mean0,checklist=checklist,check=check))
}
