write.obj <- function(obj,filename="default")
{
    obj <- format(obj,scientific=FALSE,trim=TRUE, justify = "none")
    write.table(obj,file=paste(filename,".obj",sep=""),quote=F,row.names = FALSE,col.names = FALSE,na="")
}
