anonymize <- function(data,path=NULL,dest.path=NULL,suffix=".ply",split="_",remove=NULL)
{
  if (length(dim(data)) == 3)
    {
      orig <- dimnames(data)[[3]]
    }
  else
    {
      orig <- rownames(data)
    }
  alias.tmp <- strsplit(orig,split=split)
  checklen <- unlist(lapply(alias.tmp,length))
  faclen <- checklen[1]
  checklen <- expand.grid(checklen,checklen)
  lenbool <- identical(checklen[,1],checklen[,2])
  if (!lenbool)
    {stop("please check dimnames for similar entrie numbers")
   }
  
  alias.frame <- data.frame(matrix(NA,length(orig),faclen))
  for (i in 1:length(orig)){alias.frame[i,] <- alias.tmp[[i]]}
  alias <- sprintf("%04d",1:length(orig))
  alias.frame[,remove[1]] <- alias
  if (length(remove) > 1)
    {
      alias.frame <- alias.frame[,-(remove[-1])]
    }
  alias <- NULL; for (i in 1:length(orig)){alias[i] <- paste(alias.frame[i,],collapse="_")}
  anonymkey <- data.frame(orig,alias)
  if (!is.null(path))
  {
    for ( i in 1:length(orig))
      {
        system(paste("cp ",path,"/",orig[i],suffix," ", dest.path,"/",alias[i],suffix,sep=""))
      }
  }
  if (length(dim(data)) == 3)
    {
      dimnames(data)[[3]] <- alias
    }
  else
    {
      rownames(data) <- alias
    }
  
  
  out <- list(data=data,anonymkey=anonymkey)
  return(out)
}

