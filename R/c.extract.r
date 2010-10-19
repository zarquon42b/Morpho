c.extract<-function(pts.file)
{	x<-pts.file	
	allnames<-row.names(x)
	cs<-grep("C",allnames)
	cnames<-row.names(x)[cs]	
	t<-levels(as.factor(substr(cnames,1,4)))
	tl<-length(t)
	out<-list()	
		for (i in 1:tl)
		{out[[t[i]]]<-grep(t[i],allnames)}
	
	return(out)

}
