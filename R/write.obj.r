write.obj <- function(obj,filename="default")
{write.table(obj,file=paste(filename,".obj",sep=""),quote=F,row.names = FALSE,col.names = FALSE,na="")}