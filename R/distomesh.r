distomesh<-function(dataarray,meshlist,landm,margin,reflect=TRUE,tol=5)
{

      if (is.data.frame(meshlist)==T)
      {meshlist<-list(meshlist)}

      k<-dim(dataarray)[1]
      m<-dim(dataarray)[2]
      #tol<-tol^2
      if (length(dim(dataarray))<3 )
        {
          dataarray<-array(dim=c(k,m,1),dataarray)
        }
      
      if (is.list(landm)==F )
        {
          landm<-list(landm)
        }
      n<-dim(dataarray)[3]
      nm<-length(meshlist)
      euklid.array<-array(dim=c(k,nm,n),NA)
      out<-matrix(NA,n,nm)
      test.vector<-c(rep(NA,n))
      fit.list<-list(numeric(0))
      sur.list<-list(numeric(0))
      land.list<-list(numeric(0))
      for (j in 1:nm)
          {
          #proc.test<-rotonto(landm[[j]],dataarray[,,1],scale=F,reflect=reflect)
          #if ( max(abs(proc.test$Ahat-landm[[j]])) > 1e-12)
          #{meshlist[[j]]<-warp_obj(meshlist[[j]],landm[[j]],proc.test$Ahat)
          #  cat(paste("mesh rotation for mesh #", j ,"executed","\n"))}
          
          vert<-as.matrix(subset(meshlist[[j]],meshlist[[j]][,1]=="v")[,2:4])
          nv<-dim(vert)[1]
          lmlist<-list(numeric(0))
          fit.list[[j]]<-dataarray
          sur.list[[j]]<-dataarray
          
          for (i1 in 1:k)
            {  lmlist[[i1]]<-numeric(0)
               a1<-(abs(vert[,1]-landm[[j]][i1,1]))<tol
               a2<-(abs(vert[,2]-landm[[j]][i1,2]))<tol
               a3<-(abs(vert[,3]-landm[[j]][i1,3]))<tol
               lmlist[[i1]]<-which((a1*a2*a3) == 1)
                
               #print(length(lmlist[[i1]]));print(dim(vert))              #for (j1 in 1:nv)
              #{  p<-length(lmlist[[i1]])
              #
              #    if (sum((vert[j1,]-landm[[j]][i1,])^2) < tol)
             #     {lmlist[[i1]][p+1]<-j1}
             # }
              if (length(lmlist[[i1]])==0)
                {
                  
                  cat("Landmark #",i1,"further from next point on surface than tolerated. please expand tolerance..\n")
                  stop("at least one Landmark further from next point on surface than tolerated. please expand tolerance to an adequate threshold.")
                }
            }
            
          
          for (i in 1:n)
          { cat(i)
          #Procrustes Fit ohne Skalierung der zu vergleichenden Matrizen
          proc.a<-rotonto(landm[[j]],dataarray[,,i],scaling=F,signref=F)
          fit.list[[j]][,,i]<-proc.a$yrot
          land.list[[j]]<-landm[[j]]
          testa<-proc.a$X

          
             for (t2 in 1:k)
                {
                #dif<-apply(vert[lmlist[[t2]],],1,function(x){x-proc.a$yrot[t2,]})
                #dif1<-diag(crossprod(dif))
                  
                 # spheres3d(vert[lmlist[[t2]],],col=3,radius=0.5)
                
               # m1<-which(dif1==min(dif1))
                
                for (t1 in 1: length(lmlist[[t2]]))
                 { 
                    if (sum((vert[lmlist[[t2]][t1],]-proc.a$yrot[t2,])^2) <= sum((testa[t2,]-proc.a$yrot[t2,])^2))
                    {testa[t2,]<-vert[lmlist[[t2]][t1],]}
                 }
                #testa[t2,]<-vert[lmlist[[t2]][m1],] 
                }
          sur.list[[j]][,,i]<-testa
          dif.proc<-testa-proc.a$yrot



                          
          # Berechnung der Euklid Distanzen zwischen den analogen Landmarks und Eintrag in array
          euklid.array[,j,i]<- sqrt(diag(tcrossprod(dif.proc)))
          # Abgleich der Distanzen mit dem Vorgegebenen Wert

      if (max(euklid.array[,j,i])<=margin)
      {out[i,j]<-0}
      else
      {out[i,j]<-1}

   }        
        # Berechnung des Anteils der Individuen mit Abweichungen <= margin
        
          
        
  }
        test.vector<-apply(out,1,prod)
        
        accuracy<-matrix(round(100-(sum(test.vector)/n*100),digits=2),1,1)
        colnames(accuracy)<-"%"
        rownames(accuracy)<-"accuracy"

  return(list(dif=dif,out=out,euklid.array=euklid.array,test.vector=test.vector,accuracy=accuracy,fit.list=fit.list,land.list=land.list,sur.list=sur.list,rot.mesh=meshlist))
}                                                                                              