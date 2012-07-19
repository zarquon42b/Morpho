regdist.raw<-function(dataarray,plot=TRUE,main="",rho="angle",dist.mat.out=FALSE)
{     proc<-procSym(dataarray)
      x<-proc$rotated
      n<-dim(x)[3]
      m<-dim(x)[2]
      k<-dim(x)[1]
      y<-proc$orpdata





      qm<-dist(t(matrix(x,k*m,n)))  #calc  dist. between rotated config
      procdis<-sum(qm^2)/n


        procdistmat<-matrix(NA,n,n) #calc rho from angle between rotated configs
          for (i in 1:n)
           {for (j in 1:n)
           if (rho=="riemdist"){{procdistmat[i,j]<-riemdist(x[,,i],x[,,j])}}  # riemann dist.
           else if (rho=="angle"){{procdistmat[i,j]<-angle.calc(x[,,i],x[,,j])$rho}}
            }

           if(rho=="sindist")
            {procvec<-asin(qm)}

           else
            {
              procvec<-as.dist(procdistmat)
            }

           procdis2<-sum(procvec^2)/n

      em<-dist(t(matrix(y,k*m,n)))
      euvec<-(em)
      eudis<-sum(euvec^2)/n



      correlation<-cor(euvec,procvec)^2
      if(plot==TRUE){
      plot(euvec,procvec,asp=1,xlab="euclid. dist. in tangentspace",ylab=paste("rho as",rho),main=main)
      abline(0,1,col="grey50")
      }
      if (dist.mat.out==TRUE)
        {
          return(list(cor=correlation,procSS=procdis,tanSS=eudis,rhoSS=procdis2,euc.dist=em,proc.dist=procvec))
        }
      
      else
        {
          return(list(cor=correlation,procSS=procdis,tanSS=eudis,rhoSS=procdis2))
         }
}
