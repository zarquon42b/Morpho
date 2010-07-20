read.csv.folder<-function(folder,x,y,rownames,header=TRUE,dec=".",sep=";",pattern="csv")
{	
	file.ext<-paste(".",pattern,sep="")
	name<-list.files(folder,pattern=file.ext)
	
		
	ln<-length(name)
	arr<-array(NA,dim=c(length(x),3,ln))
	#print(dim(arr))
	
	if (is.character(x))
	for ( i in 1:ln)	
		{data<-read.table(paste(folder,name[i],sep=""),header=header,dec=dec,sep=sep)
		#print(dim(data[x,y]))
		dat<-NULL
		count<-1
		rn<-data[,rownames]
		for (j in 1:length(x))
			{
			dat[count]<-which(rn==x[j])
			count<-count+1
			}
		print(dat)
		
		 		
		#dat<-which(data[,rownames] %in% x)
		#print(dat)
		arr[,,i]<-as.matrix(data[dat,y])
		if (i ==1)
		rown<-data[dat,rownames]
		
		
		}
	else
		{
			for ( i in 1:ln)
			{data<-read.table(paste(folder,name[i],sep=""),header=header,dec=dec,sep=sep)
			#print(dim(data[x,y]))
			arr[,,i]<-as.matrix(data[x,y])
			if (i ==1)
			rown<-data[x,rownames]
			}
		}
	xlen<-length(x);ylen<-length(y)
	nas<-which(is.na(arr))
	nas<-as.integer(nas/(xlen*ylen))+1
	nas<-nas[-(which(duplicated(nas)))]
	
	
	dimnames(arr)<-list(rown,c("X","Y","Z"),sub(file.ext,"",name))
	return(list(arr=arr,NAs=nas))
	

}
