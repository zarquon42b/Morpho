#' batch import data from files
#' 
#' imports all data files contained in a specified folder.
#' 
#' 
#' @param folder character: path to folder
#' @param x either a vector specifiing which rows are to be imported, or
#' character vector containing variable names to be sought for.
#' @param y a vector specifiing, which columns of the speradsheet ist to be
#' imported.
#' @param rownames integer: specifies columns, where variable names are stored.
#' @param header logical : if spreadsheet contains header-row.
#' @param dec character: defines decimal sepearator.
#' @param sep character: defines column seperator.
#' @param pattern character: specify file format (e.g. csv).
#' @param addSpec character: add a custom specifier to the dimnames of the
#' array.
#' @param back logical: where to place the specifier.
#' @return
#' \item{arr }{array containing imported data}
#' \item{NAs }{vector containing position of observations with NAs}
#' \item{NA.list }{list: containing vectors containing information which
#' LMs are missing in which observation}
#' @author Stefan Schlager
#' @seealso \code{\link{read.table}}
#' @export
read.csv.folder <- function(folder,x,y=2:4,rownames=NULL,header=TRUE,dec=".",sep=";",pattern="csv",addSpec=NULL,back=TRUE)
{	
	if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/")
		{folder <- paste(folder,"/",sep="")
		}
	
	file.ext <- paste(".",pattern,sep="")
	name <- list.files(folder,pattern=file.ext)
	xlen <- length(x)
	ylen <- length(y)
	NA.list <- NULL	
	
	ln <- length(name)
	arr <- array(NA,dim=c(xlen,ylen,ln))
	if (is.factor(x))
		{x <- as.character(x)
		}	
	if (is.character(x)) ### check if selection contains variable names 
		for ( i in 1:ln)	
		{data <- read.table(paste(folder,name[i],sep=""),header=header,dec=dec,sep=sep)
		dat <- NULL
		count <- 1
		if (is.null(rownames))
			{stop("please specify column containing Landmark names!")
			}
		rn <- data[,rownames]
		for (j in 1:length(x))
			{check <- which(rn==x[j])
			
			if (length(check)==0)
				{warning(paste("dataset",i,"misses entry for Landmark",j))
				data[9999,y] <- rep(NA,ylen)
				dat[count] <- 9999
				
				}
			if (length(check) > 1)
				{warning(paste("dataset",i,"contains landmark #",x[j],"with the same name - first match was used."))
				dat[count] <- check[1]
				}
			else
				{empty <- which(rn==x[j])
                                 if (length(empty) !=0)
                                   {
                                     dat[count] <- which(rn==x[j])
                                   }
				}
			count <- count+1
			}
		arr[,,i] <- as.matrix(data[dat,y])
		if (i ==1)
		rown <- x
		
		}
	
	else
		{for (i in 1:ln)
			{data <- read.table(paste(folder,name[i],sep=""),header=header,dec=dec,sep=sep)
			arr[,,i] <- as.matrix(data[x,y])
			if (i ==1)
			if (is.null(rownames))
				{rown <- c(1:xlen)
				}
			else
				{rown <- data[x,rownames]
				}
			}
		}
	
	
	nas0 <- which(is.na(arr))	### check for NAs and store information about missing Landmark and individual
	nas1 <- as.integer(nas0/(xlen*ylen))+1
	nas <- nas1[-(which(duplicated(nas1)))]
        
	if (length(nas) > 0)
		{NA.list <- list()
			for (i in 1:length(nas))
			{nas2 <- nas0[which(nas1==nas[i])]%%(xlen*ylen)
			nas2 <- nas2%%xlen
			nas2 <- nas2[-which(duplicated(nas2))]
			if (0 %in% nas2)
				{nas2[which(nas2==0)] <- xlen}
                         if (length(nas2) > 0)
                           {
                             NA.list[[as.character(nas[i])]] <- sort(nas2)
                           }
                         else
                           {
                             NA.list[[as.character(nas[i])]] <- NULL
                             nas <- nas[-i]
                           }
                       }
               }
	else 
		{nas <- NULL
		}	
	
	dim3 <- NULL
	if (back)
		{dim3 <- paste(sub(file.ext,"",name),addSpec,sep="")
		}
	else
		{dim3 <- paste(addSpec,sub(file.ext,"",name),sep="")
		}
	if (ylen==2)
		{dimnames(arr) <- list(rown,c("X","Y"),dim3)
		}
	else
	
	{dimnames(arr) <- list(rown,c("X","Y","Z"),dim3)
	}
		
	return(list(arr=arr,NAs=nas,NA.list=NA.list))
	

}
