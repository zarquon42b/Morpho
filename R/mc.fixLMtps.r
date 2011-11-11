mc.fixLMtps<-function(data,comp=3,weight=TRUE)
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
		
		### rotate incomplete data onto mean ###
		lmsdat<-data
		
		#print(data.c)
		proc.c<-mc.procGPA(data.c)
		mean0<-proc.c$mshape
		
		for (i in 1:length(check))
			{
			miss<-checklist[[check[i]]]
			if (weight) ### calculate weighted estimates of missing data ###
				{
			### rotate incomplete data onto mean ###
				rotmiss<-rotonto(mean0[-miss,],data[-miss,,check[i]],scale=TRUE)$yrot
				allrot<-abind(rotmiss,proc.c$rotated[-miss,,])
			### calculate weights according to procrustes distance ###			
				wcalc<-proc.weight(allrot,comp,1,report=FALSE)
				lms<-proc.c$rotated[,,wcalc$data$nr-1]
				lm.est<-matrix(0,dim(data)[1],m)
			
			for (j in 1:comp)
				{lm.est<-lm.est+lms[,,j]*wcalc$data$weight[j]
				}
			
			tpsout<-tps3d(lm.est,lm.est[-miss,],data[-miss,,check[i]])
			}
			else
				{tpsout<-tps3d(mean0,mean0[-miss,],data[-miss,,check[i]])
				}
			#print(tpsout)
			lmsdat[,,check[i]]<-lm.est
			
			out[,,check[i]]<-tpsout
			}
		return(list(out=out,mshape=mean0,checklist=checklist,check=check,lmsdat=lmsdat))
}
