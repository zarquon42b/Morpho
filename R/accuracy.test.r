accuracy.test<-function(euklid.array,margin)
{   n<-dim(euklid.array)[3]
    m<-dim(euklid.array)[2]
    k<-dim(euklid.array)[1]
        out<-matrix(NA,n,m)
        test.vector<-c(rep(NA,n))
    analysis<-c(rep(0,k))

  for (i in 1:n)
  {for (j in 1:m)
  { if (max(euklid.array[,j,i])<margin)
  {out[i,j]<-0}
  else
  {
    out[i,j]<-1}
    for (i1 in 1:k)
      {
      if (min(euklid.array[i1,,i])>=margin)
        analysis[i1]<-analysis[i1]+1
      }
  }
  # Berechnung des Anteils der Individuen mit Abweichungen <= margin
        test.vector[i]<-prod(out[i,])
  }
        accuracy<-matrix(round(100-(sum(test.vector)/n*100),digits=2),1,1)
        colnames(accuracy)<-"%"
        rownames(accuracy)<-"accuracy"
        analysis<-as.data.frame(analysis)
        colnames(analysis)<-">= margin"
            
return(list(out=out,accuracy=accuracy,analysis=analysis))}
