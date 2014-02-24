#' Graphical interface to find outliers and/or to switch mislabeld landmarks
#' 
#' Graphical interface to find outliers and/or to switch mislabeld landmarks
#' 
#' This function performs a procrustes fit and sorts all specimen according to
#' their distances (either Procrustes or Mahalanobis-distance) to the sample's
#' consensus. It provides visual help for rearranging landmarks and/or
#' excluding outliers.
#' 
#' @param A Input k x m x n real array, where k is the number of points, m is
#' the number of dimensions, and n is the sample size.
#' @param color color of Landmarks points to be plotted
#' @param lwd linewidth visualizing distances of the individual landmarks from
#' mean.
#' @param lcol color of lines visualizing distances of the individual landmarks
#' from mean.
#' @param mahalanobis logical: use mahalanobis distance to find outliers.
#' @param PCuse integer: Restrict mahalanobis distance to the first n Principal
#' components.
#' @param text logical: if \code{TRUE}, landmark labels (rownumbers) are
#' displayed
#' @return
#' \item{data.cleaned }{array (in original coordinate system) containing
#' the changes applied and outliers eliminated}
#' \item{outlier }{vector with integers indicating the positions in the
#' original array that have been marked as outliers}
#' \item{dist.sort }{table showing the distance to mean for each
#' observation - decreasing by distance}
#' \item{type }{what kind of distance was used}
#' @author Stefan Schlager
#' @seealso \code{\link{typprob}},\code{\link{typprobClass}}
#' 
#' @examples
#' 
#' data(boneData)
#' ## look for outliers using the mahalanobis distance based on the first
#' # 10 PCscores
#' # to perform the example below, you need,of course, uncomment the answers
#' \dontrun{
#' outliers <- find.outliers(boneLM, mahalanobis= TRUE, PCuse=10)
#' # n # everything is fine 
#' # n # proceed to next
#' # s # let's switch some landmarks (3 and 4)
#' # 3
#' # 4
#' # n # we are done
#' # y # yes, because now it is an outlier
#' # s #enough for now
#' }
#' 
#' @export
find.outliers <- function(A,color=4,lwd=1,lcol=2,mahalanobis=FALSE,PCuse=NULL, text=TRUE)
{   	
    raw <- A
    n <- dim(A)[3]
    k <- dim(A)[1]
    m <- dim(A)[2]
    A <- ProcGPA(A)
    disType <- "Procrustes"
###from here on the same as find.outliers ###		
    rho <- NULL
    if (!mahalanobis) {
        a.list <- as.list(1:n)
        rhofun <- function(i)
            {
                rho <- angle.calc(A$rotated[,,i],A$mshape)
                return(rho)          	
            }
        rho <- unlist(lapply(a.list,rhofun))
        A$rho <- rho
    } else {
        cat("mahalanobis distance is used\n")
        PCs <- prcomp(vecx(A$rotated))$x
        if (!is.null(PCuse))
            PCs <- PCs[,1:PCuse]
        A$rho <- mahalanobis(PCs,cov=cov(PCs),center=0)
        disType <- "Mahalanobis D^2"     
    }

    if (is.null(dimnames(raw)[[3]]))
        dimnames(raw)[[3]] <- rep("",n)
    
    disti <- data.frame(c(1:n),A$rho,dimnames(raw)[[3]])
    disti.sort <- disti[order(disti[,2],decreasing=T),]
    colnames(disti.sort)[1] <- "# in array"
    rownames(disti.sort) <- c(1:n)
    outlier <- NULL
    t1 <- 1
    lo <- 1
    if (m==3)
        op <- open3d()	
    
    while (t1 <= n ) {
        if (m==3)
            .difplotLM(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE, text=text)
        else 
            .difplotLM2D(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,main=disti.sort[t1,1], text=text)

        cat(paste("outlier #",t1,": ",disti.sort[t1,1]," - ",disti.sort[t1,3],"     ",disType," dist. to mean: ",disti.sort[t1,2],"\n",sep=""))
        
        if (disti.sort[t1,1] %in% outlier) {
            answer <- substr(readline(" already added to outlierlist! remove from list (y/N/s)?\ny=yes,n=no,s=switch landmarks: "), 1L,1L)
           
            while (!(answer %in% c("y","Y","s","S","N","n"))) {
                answer <- substr(readline(" already added to outlierlist! remove from list (y/N/s)?\ny=yes,n=no,s=switch landmarks: "), 1L,1L)
            }
### keep in outlier list		
            if (answer == "n" || answer == "N"){
                lo <- lo+1
                ## switching landmark routine starting ###	 	
            } else if (answer == "s" || answer == "S") {
                loop0="y"
                while(loop0=="y" || loop0 =="Y") {
                    answer1a <- as.integer(readline("select  two Landmarks to switch position\ninsert first: "))
                    while (is.na(answer1a)|| answer1a > k || answer1a < 1) {
                        answer1a <- as.integer(readline(paste("please enter integer <",k,"!\ninsert first: ")))
                    }
                    answer1b <- as.integer(readline("insert second: "))
                    while (is.na(answer1b) || answer1b > k || answer1b < 1) {
                        answer1b <- as.integer(readline(paste("please enter integer <",k,"!\ninsert second: ")))
                    }
                    
### switch rows of selected landmarks ##
                    raw[c(answer1a,answer1b),,disti.sort[t1,1]] <- raw[c(answer1b,answer1a),,disti.sort[t1,1]]
                    A$rotated[c(answer1a,answer1b),,disti.sort[t1,1]] <- A$rotated[c(answer1b,answer1a),,disti.sort[t1,1]]
                    rho.new <- angle.calc(A$rotated[,,disti.sort[t1,1]],A$mshape)
                    
                    if (m==3) {
                        rgl.clear()
                        .difplotLM(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE, text=text)
                    } else {
                        .difplotLM2D(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,main=disti.sort[t1,1],text=text)
                    }
                    cat(paste("new distance to mean:",rho.new,"\n"))
                    loop0 <- substr(readline("switch more (y/N)? "),1L,1L)
                    while(loop0 != "y" && loop0 != "Y" && loop0 != "n" && loop0 != "N" ) {
                        loop0 <- substr(readline("yes or no? "),1L,1L)
                    }
                }
                answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
                while (!(answer %in% c("y","Y","N","n")))
                    {answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
                 }   
                
### remove from outlier list ###
            } else if (answer == "y" || answer == "Y") {
                outlier <- outlier[-which(outlier==disti.sort[t1,1])]
                cat(paste("outlier #",t1,": ",disti.sort[t1,1],"removed from outlierlist.\n"))
                lo <- lo+1
            }	
        } else {
            answer <- substr(readline(" add to outlierlist (y/N/s)?\ny=yes,n=no,s=switch landmarks: "), 1L,1L)
            while (!(answer %in% c("y","Y","s","S","N","n"))) {
                answer <- substr(readline(" add to outlierlist (y/N/s)?\ny=yes,n=no,s=switch landmarks: "), 1L,1L)
            }
### add to outlier list		
            if (answer == "y" || answer == "Y") {
                outlier[lo] <- disti.sort[t1,1]
                lo <- lo+1
                
### switching landmark routine starting ###	 	
            } else if (answer == "s" || answer == "S") {
                loop0="y"
                while(loop0=="y" || loop0 =="Y") {
                    answer1a <- as.integer(readline("select  two Landmarks to switch position\ninsert first: "))
                    while (is.na(answer1a)|| answer1a > k || answer1a < 1) {
                        answer1a <- as.integer(readline(paste("please enter integer <",k,"!\ninsert first: ")))
                    }
                    answer1b <- as.integer(readline("insert second: "))
                    while (is.na(answer1b) || answer1b > k || answer1b < 1) {
                        answer1b <- as.integer(readline(paste("please enter integer <",k,"!\ninsert second: ")))
                    }
                    
### switch rows of selected landmarks ##
                    raw[c(answer1a,answer1b),,disti.sort[t1,1]] <- raw[c(answer1b,answer1a),,disti.sort[t1,1]]
                    A$rotated[c(answer1a,answer1b),,disti.sort[t1,1]] <- A$rotated[c(answer1b,answer1a),,disti.sort[t1,1]]
                    rho.new <- angle.calc(A$rotated[,,disti.sort[t1,1]],A$mshape)
                    if (m==3) {
                        rgl.clear()
                        .difplotLM(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,rgl.new=FALSE, text=text)
                    } else {
                        .difplotLM2D(A$mshape,A$rotated[,,disti.sort[t1,1]],color=color,lwd=1,lcol=lcol,main=disti.sort[t1,1], text=text)
                    }
                    cat(paste("new distance to mean:",rho.new,"\n"))
                    loop0 <- substr(readline("switch more (y/N)? "),1L,1L)
                    while(loop0 != "y" && loop0 != "Y" && loop0 != "n" && loop0 != "N" ) {
                        loop0 <- substr(readline("yes or no? "),1L,1L)
                    }
                }
                answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
                while (!(answer %in% c("y","Y","N","n"))) {
                    answer <- substr(readline("do you still want to add to outliers (y/N)? "), 1L,1L)
                 }
                if (answer == "y" || answer == "Y") {
                    outlier[lo] <- disti.sort[t1,1]
                    lo <- lo+1
                }	
            }
        }	
### go on or halt ###
        answer1 <- substr(readline(" next/previous/stop (n/p/s)?  "), 1L,1L)
        if (!(answer1 %in% c("P","p","S","s"))) {
            if (m==3) {
                rgl.clear()
                rgl.bg(color = "white")
            }
            t1 <- t1+1          
        } else if (answer1 %in% c("P","p")) {
            if (m==3) {
                rgl.clear()
                rgl.bg(color = "white")
            }
            if (t1 !=1) {
                t1 <- t1-1
            } else 
                cat("already at top of the line!\n")
        } else if (answer1 == "s" || answer1 == "S")
            break 
    }
### remove outlier from array
    if (! is.null(outlier))
        raw <- raw[,,-outlier]
    invisible(list(data.cleaned=raw,outlier=outlier,dist.sort=disti.sort,type=disType))
}
