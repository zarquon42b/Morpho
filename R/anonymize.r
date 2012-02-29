anonymize <- function(data,path=NULL,dest.path=NULL,ext=".ply",split="_",remove=NULL,levels=TRUE,prefix=NULL,suffix=NULL,sample=TRUE)
{
  if (length(dim(data)) == 3)
    {
      orig <- dimnames(data)[[3]]
    }
  else if (length(dim(data)) == 2)
    {
      orig <- rownames(data)
    }
  else

    {orig <- data
   }
      
  alias.tmp <- strsplit(orig,split=split)
  checklen <- unlist(lapply(alias.tmp,length))
  faclen <- checklen[1]
  checklen <- expand.grid(checklen,checklen)
  lenbool <- identical(checklen[,1],checklen[,2])
  if (!lenbool)
    {
      stop("please check dimnames for similar entry numbers")
    }
  
  alias.frame <- data.frame(matrix(NA,length(orig),faclen))
  
  for (i in 1:length(orig)){alias.frame[i,] <- alias.tmp[[i]]}
  
  orig.frame <- alias.frame
  if (length(remove) > 1)
    {
      facs <- apply(orig.frame[,remove],1,paste,collapse="_")
    }
  else
    { facs <- orig.frame[,remove]
    }
  if (levels)
    {
      levorig <- levels(as.factor(facs))
     
        
      alias.lev <- paste(prefix,sprintf("%04d",1:length(levorig)),suffix,sep="")
       if (sample)
         {
           alias.lev <- sample(alias.lev)
         }
      alias <- rep(0,length(orig))
      for (i in 1:length(levorig))
        {
          alias[facs==levorig[i]] <- alias.lev[i]
        }
    }
  else
    {
      alias <- sample(sprintf("%04d",1:length(orig)))
    }
  
  alias.frame[,remove[1]] <- alias
  if (length(remove) > 1)
    {
      alias.frame <- alias.frame[,-(remove[-1])]
    }

  alias <- NULL
  if (!is.null(dim(alias.frame)))
    {
      for (i in 1:length(orig)){alias[i] <- paste(alias.frame[i,],collapse="_")}
    }
  else
    {
      for (i in 1:length(orig))
        {
          alias[i] <-alias.frame[i]
        }
    }
  anonymkey <- data.frame(orig,alias)
  if (!is.null(path))
    {
      for ( i in 1:length(orig))
        {
          system(paste("cp ",path,"/",orig[i],ext," ", dest.path,"/",alias[i],ext,sep=""))
        }
    }
  if (length(dim(data)) == 3)
    {
      dimnames(data)[[3]] <- alias
    }
  else if  (length(dim(data)) == 2)
    {
      rownames(data) <- alias
    }
  else
    {
      data <- alias
    }
  
  out <- list(data=data,anonymkey=anonymkey)
  return(out)
}
