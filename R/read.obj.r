read.obj<-function(file)
{   
	test<-read.table(file,sep=" ",nrows=50)
	if(length(grep(",",test[30:50,2])) != 0)
		{out<-read.table(file,sep=" ",dec=",")}
	
	else {out<-read.table(file,sep=" ")}
   	#vn<-(grep("vn",out[,1]))
    	#	if (length(vn)!=0)
      	#	{out<-out[-vn,]}
    
    out
}
