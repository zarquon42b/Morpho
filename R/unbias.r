unbias<-function(good,wrong,refl,refr,left,right,directional=TRUE)
{		
		
		vectl<-good[refl,]-wrong[refl,]
		if (length (refl)>1)
			{vectl<-apply(vectl,2,mean)
			}
		print(vectl)
		addl<-matrix(vectl,dim(good[left,])[1],dim(good[left,])[2],byrow=T)
		vectr<-good[refr,]-wrong[refr,]
		if (length (refr)>1)
			{vectr<-apply(vectr,2,mean)
			}
		addr<-matrix(vectr,dim(good[right,])[1],dim(good[right,])[2],byrow=T)
		out<-wrong
		out[left,]<-wrong[left,]+addl
		out[right,]<-wrong[right,]+addr
	return(out)
}
		
