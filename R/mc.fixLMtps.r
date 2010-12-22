mc.fixLMtps<-function(data)
{		n<-dim(data)[3]
		k<-dim(data)[1]
		m<-dim(data)[2]
		checklist<-list()
		checkvec<-rep(0,n)
		out<-data
	### check for missing landmarks ###
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
		#### calc mean of complete configs ###
		check<-which(checkvec==1)
		data.c<-data[,,-check]
		#print(data.c)
		mean0<-mc.procGPA(data.c)$mshape
		
		for (i in 1:length(check))
			{miss<-checklist[[check[i]]]
			out[,,check[i]]<-mc.tps3d(mean0,mean0[-miss,],data[-miss,,check[i]])
			}
		return(list(out=out,mshape=mean0,checklist=checklist,check=check))
}
