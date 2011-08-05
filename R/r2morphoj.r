r2morphoj <- function(x,file,id.string=NULL)
  {
    x <- vecx(x,byrow=TRUE)
    if (is.null(id.string))
      {
        id.string <- dimnames(x)[[3]]
      }
    out <- data.frame(id.string,x,row.names=id.string)
    write.table(out,file=file,quote=F,row.names=FALSE,sep="\t")
    
    
  }
