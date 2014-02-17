#' Export data to MorphoJ and Morphologika
#' 
#' Export data to MorphoJ and Morphologika
#' 
#' 
#' @title Export data to MorphoJ and Morphologika
#' @param x 3-dimensionla array containing landmark data. E.g. the input/output
#' from \code{\link{procSym}}.
#' @param file character: name the output file
#' @param id.string a string with ids or factors to append
#' @param labels character vector specify labels to create for Morphologika
#' @param labelname character: name the labels for Morphologika.
#' @param \dots unused at the moment
#' 
#' @examples
#' 
#' library(shapes)
#' r2morphoj(gorf.dat,file="gorf.dat")
#' 
#' data <- bindArr(gorf.dat, gorm.dat, along=3)
#' datalabels <- c(rep("female",dim(gorf.dat)[3]),
#' rep("male",dim(gorm.dat)[3]))
#' labelname <- "sex"
#' r2morphologika(data, labels=datalabels, labelname= labelname, file="data.dat")
#' 
#' @rdname r2morphoj
#' @export
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
