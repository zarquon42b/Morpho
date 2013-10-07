read.obj <- function(file)
{   
	test <- read.table(file, nrows=50)
	if(length(grep(",",test[30:50,2])) != 0)
            out <- read.table(file,dec=",")
        else
            out <- read.table(file)
   	#vn <- (grep("vn",out[,1]))
    	#	if (length(vn)!=0)
      	#	{out <- out[-vn,]}
    
    out
}
