#' @rdname r2morphoj
#' @export
r2morphologika <- function(x,file=file,labels=NULL,labelname=NULL,...)
  {
    n <- dim(x)[3]
    m <- dim(x)[2]
    k <- dim(x)[1]
    idnames <- dimnames(x)[[3]]
    if (is.null(idnames))
      idnames <- paste("specimen",1:n)

    ## start writing to file
    cat("[individuals]\r\n",file = file,sep="")
    cat(paste(n,"\r\n",sep=""),file = file,append=TRUE,sep="")
    cat("[landmarks]\r\n",file = file,append=TRUE,sep="")
    cat(paste(k,"\r\n",sep=""),file = file,append=TRUE,sep="")
    cat("[dimensions]\r\n",file = file,append=TRUE,sep="")
    cat(paste(m,"\r\n",sep=""),file = file,append=TRUE,sep="")
    cat("[names]\r\n",file = file,append=TRUE,sep="")
    cat(paste(idnames,"\r\n",sep=""),file = file,append=TRUE,sep="")
    if (!is.null(labels))
      {
        cat("[labels]\r\n",file = file,append=TRUE,sep="")
        if (is.null(labelname))
          labelname <- as.character(bquote(labels))
        cat(paste(labelname,"\r\n",sep=""),file = file,append=TRUE,sep="")
        cat("[labelvalues]\r\n",file = file,append=TRUE,sep="")
        cat(paste(labels,"\r\n",sep=""),file = file,append=TRUE,sep="")
      }
    
    cat("[rawpoints]\r\n",file = file,append=TRUE,sep="")
    for (i in 1:n)
      {
        cat(paste("\'#",i,"\r\n",sep=""),file = file,append=TRUE,sep="")
        write.table(format(x[,,i], scientific = F, trim = T), file = file, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, na = "", eol="\r\n")
        
     }
        
  }

     
