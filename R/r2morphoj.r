r2morphoj <- function(x,file,id.string=NULL)
  {
    x <- vecx(x,byrow=TRUE)
    if (is.null(id.string))
      {
          if (is.null(rownames(x)))
              id.string <- paste(1:nrow(x))
          else
              id.string <- rownames(x)
      }
    out <- data.frame(id.string,x,row.names=NULL)
    write.table(out,file=file,quote=F,row.names=FALSE,sep="\t")
}
