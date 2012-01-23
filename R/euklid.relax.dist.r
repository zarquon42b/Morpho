euklid.relax.dist<-function(DataArray,RefArray,margin,SMvector=0,outlines=0,deselect=F)
{   library(shapes)
    
    n<-dim(DataArray)[3]
    k<-dim(DataArray)[1]
    m<-dim(DataArray)[2]     
        
    if (length(dim(RefArray))<3 )
    {
      RefArray<-array(dim=c(k,m,1),RefArray)
    }  
    
    l<-dim(RefArray)[3]
    euklid.frame<-array(dim=c(k,l,n),NA)
    out<-matrix(NA,n,l)
    test.vector<-c(rep(NA,n))
    dataslide<-DataArray    
    fit.array<-DataArray   
  for (i in 1:n)
    
  {    
        
    for (j in 1:l)
          {  
        if (SMvector[1]!=0)
          {
            if (m==3)
              {L<-CreateL(RefArray[,,j])}
            else 
              {L<-CreateL2D(RefArray[,,j])} 
              #for (j in 1:n)
             
              U<-calcTang_U(DataArray[,,i],SMvector=SMvector,outlines=outlines,deselect=deselect)
              dataslide[,,i]<-calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m)$Gamatrix
            
          }
          #Procrustes Fit ohne Skalierung der zu vergleichenden Matrizen 
          proc.a<-procOPA(RefArray[,,j],dataslide[,,i],scale=F,reflect=T)
          fit.array[,,i]<-proc.a$Bhat
          dif.proc<-proc.a$Ahat-proc.a$Bhat
          

          # Berechnung der Euklid Distanzen zwischen den analogen Landmarks und Eintrag in array
          euklid.frame[,j,i]<- sqrt(diag(tcrossprod(dif.proc)))   
          
  
      # Abgleich der Distanzen mit dem Vorgegebenen Wert
        
      if (max(euklid.frame[,j,i])<=margin)
      {out[i,j]<-0}
      else
      {out[i,j]<-1}
       
     }
        # Berechnung des Anteils der Individuen mit Abweichungen <= margin
        test.vector[i]<-prod(out[i,])
  }
        accuracy<-matrix(round(100-(sum(test.vector)/n*100),digits=2),1,1)
        colnames(accuracy)<-"%"
        rownames(accuracy)<-"accuracy"

  return(list(out=out,euklid.array=euklid.frame,test.vector=test.vector,accuracy=accuracy,fit.array=fit.array))
}