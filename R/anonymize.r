#' Replace ID-strings of data and associated files.
#' 
#' Replace ID-strings with for digits - e.g. for blind observer error testing.
#' 
#' 
#' @param data Named array, matrix or vector containing data.
#' 
#' @param remove integer: which entry (separated by \code{split}) of the name
#' is to be removed
#' @param path Path of associated files to be copied to renamed versions.
#' @param dest.path where to put renamed files.
#' @param ext file extension of files to be renamed.
#' @param split character: by which to split specimen-ID
#' @param levels logical: if a removed entry is to be treated as a factor. E.g.
#' if one specimen has a double entry, the anonymized versions will be named
#' accordingly.
#' @param prefix character: prefix before the alias string.
#' @param suffix character: suffix after the alias ID-string.
#' @param sample logical: whether to randomize alias ID-string.
#' @return
#' \item{data }{data with names replaced}
#' \item{anonymkey }{map of original name and replaced name}
#' 
#' @examples
#' 
#' anonymize(iris,remove=1)
#' 
#' 
#' @export
anonymize <- function(data ,remove, path=NULL,dest.path=NULL,ext=".ply",split="_",levels=TRUE,prefix=NULL,suffix=NULL,sample=TRUE)
{
    if (length(dim(data)) == 3)
        orig <- dimnames(data)[[3]]
    else if (length(dim(data)) == 2)
        orig <- rownames(data)
    else
        orig <- data
    
    if (is.null(orig))
        stop("data is not named")
    alias.tmp <- strsplit(orig,split=split)
    checklen <- unlist(lapply(alias.tmp,length))
    faclen <- checklen[1]
    checklen <- expand.grid(checklen,checklen)
    lenbool <- identical(checklen[,1],checklen[,2])
    if (!lenbool)
        stop("please check dimnames for similar entry numbers")
    
    alias.frame <- data.frame(matrix(NA,length(orig),faclen))
    
    for (i in 1:length(orig)){alias.frame[i,] <- alias.tmp[[i]]}
    
    orig.frame <- alias.frame
    if (length(remove) > 1)
        facs <- apply(orig.frame[,remove],1,paste,collapse="_")
    else
        facs <- orig.frame[,remove]
    
    if (levels)
        {
            levorig <- levels(as.factor(facs))
            alias.lev <- paste(prefix,sprintf("%04d",1:length(levorig)),suffix,sep="")
            if (sample)
                alias.lev <- sample(alias.lev)
            
            alias <- rep(0,length(orig))
            for (i in 1:length(levorig))
                alias[facs==levorig[i]] <- alias.lev[i]
        }
    else
        {
            if (sample)
                alias <- sample(sprintf("%04d",1:length(orig)))
            else
                alias <- sprintf("%04d",1:length(orig))
        }
    
    alias.frame[,remove[1]] <- alias
    if (length(remove) > 1)
        alias.frame <- alias.frame[,-(remove[-1])]
    
    alias <- NULL
    if (!is.null(dim(alias.frame)))
        for (i in 1:length(orig)){alias[i] <- paste(alias.frame[i,],collapse="_")}
    else
        {
            for (i in 1:length(orig))
                alias[i] <- alias.frame[i]
        }
    anonymkey <- data.frame(orig,alias)
    if (!is.null(path))
        {
            for ( i in 1:length(orig))
                {
                    dir.create(dest.path, showWarnings = F)
                    system(paste("cp ",path,"/",orig[i],ext," ", dest.path,"/",alias[i],ext,sep=""))
                }
        }
    if (length(dim(data)) == 3)
        dimnames(data)[[3]] <- alias
    
    else if  (length(dim(data)) == 2)
        rownames(data) <- alias
    
    else
        data <- alias
    
    
    out <- list(data=data,anonymkey=anonymkey)
    return(out)
}
