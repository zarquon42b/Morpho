proc.weight<-function(data,number,ref,report=TRUE,reg=0,log=FALSE,mahalanobis=TRUE)
{	
	col<-3
        rho<-0
        if (length(dim(data))==3)
          {
            l<-dim(data)[3]
            for(i in 1:l){rho[i]<-angle.calc(data[,,ref],data[,,i])$rho}
            if (is.null(dimnames(data)[[3]]))
              {id<-as.character(c(1:l))
             }
            else
              {id<-dimnames(data)[[3]]}
          }
        else
          {
            l<-dim(data)[1]
            if (mahalanobis)
              {
                covtmp <- cov(data)
                if (reg != 0)
                  {
                    {
                      eig <- eigen(covtmp,symmetric=TRUE)
                      covtmp <- t(eig$vectors)%*%diag(eig$values+reg)%*%eig$vectors
                    }
                  }
              }
            else
              {
                covtmp <- diag(ncol(data))
              }
            rho <- mahalanobis(data,cov=covtmp,center=data[ref,])
            

            if (is.null(dimnames(data)[[1]]))
              {id<-as.character(c(1:l))
             }
            else
              {id<-dimnames(data)[[1]]}
          }
        if (log)
          {
            rho <- log(rho)
          }
        nr<-c(1:l)
        data<-data.frame(nr,id,rho)
        dat.sort.i<-data[order(data[,col]),]
        dat.which<-dat.sort.i[2:(number+1),]
        
        all<-sum(dat.which[,col])
        share<-0
        
        for (i in 1:length(dat.which[,col]))
          {share[i]<-dat.which[i,col]/all
           
         }
        weight<-sort(share,decreasing=T)
        out<-data.frame(dat.which,weight)
        if (report)
          {
            cat(paste("  reference was",id[ref],"\n"))
          }
        return(list(data=out,reference=id[ref],rho.all=data))
        
}
