#' plot histogram for multiple groups.
#' 
#' plot a histogram for multiple groups, each group colored individually
#' 
#' Just a wrapper for the function hist from the "graphics" package
#' 
#' @param data vector containing data.
#' @param groups grouping factors
#' @param main,xlab,ylab these arguments to title have useful defaults here.
#' @param col vector containing color for each group. If NULL, the function
#' "rainbow" is called.
#' @param alpha numeric between 0 and 1. Sets the transparency of the colors
#' @param breaks one of: \itemize{
#' \item a vector giving the breakpoints between histogram cells,
#' \item a single number giving the number of cells for the histogram,
#' \item a character string naming an algorithm to compute the number of cells (see \sQuote{Details}),
#' \item a function to compute the number of cells.  } In the last three cases the number is a suggestion only.
#' @param legend logical: if TRUE, a legend is plotted
#' @param legend.x x position of the legend from the upper left corner
#' @param legend.y y position of the legend from the upper left corners
#' @param legend.pch integer: define the symbol to visualise group colors
#' (\code{\link{points}})
#' @param freq logical: if TRUE, the histogram graphic is a representation of
#' frequencies, the counts component of the result; if FALSE, probability
#' densities are plotted for each group.
#' @author Stefan Schlager
#' @seealso \code{\link{hist}}
#' 
#' @examples
#' 
#' data(iris)
#' histGroup(iris$Petal.Length,iris$Species)
#' 
#' 
#' @export
histGroup <- function(data,groups, main=paste("Histogram of" , dataname),xlab=dataname,ylab,col=NULL, alpha=0.5,breaks="Sturges",legend=TRUE,legend.x=80,legend.y=80,legend.pch=15,freq=TRUE)
  {
    out <- list()
    dataname <- paste(deparse(substitute(data), 500), collapse="\n")
    
    histo <- hist(data,plot=FALSE,breaks=breaks)
    

    if(!is.factor(groups))
      {
        groups <- as.factor(groups)
      }
    lev <- levels(groups)
    nlev <- length(lev)
    if (is.null(col))
      colo <- rainbow(nlev,alpha=alpha)
    else
      {
        if (length(col) != nlev)
          stop("length of 'col' must match number of groups")
        else if (is.character(col))
          colo <- col
        else
          {
            rgbfun <- function(x)
              {
                alpha <- alpha*255
                x <- rgb(x[1],x[2],x[3],maxColorValue = 255,alpha=alpha)
                return(x)
              }
            colo <- apply(col2rgb(col),2,rgbfun)
            
          }
      }
    testrun <- 0
    for( i in 1:nlev)
     {if(freq)
       testrun[i] <- max(hist(data[groups==lev[i]],breaks=histo$breaks,plot=F)$counts)
      else
        testrun[i] <- max(hist(data[groups==lev[i]],breaks=histo$breaks,plot=F)$density)
        
     }
    ylim <- max(testrun)
    ylim <- ylim+0.15*ylim
    out[[1]] <- hist(data[groups==lev[1]],breaks=histo$breaks,col=colo[1],main=main,xlab=xlab,ylab=ylab,ylim=c(0,ylim),freq=freq)
    for (i in 2:nlev)
      {
        out[[i]] <- hist(data[groups==lev[i]],breaks=histo$breaks,col=colo[i],add=T,freq=freq)
      }
    if (legend)
      {
        tmp <- 0
        tmp[1] <- grconvertX(legend.x, 'device')
        tmp[2] <- grconvertY(legend.y, 'device') 
      legend(tmp[1],tmp[2],pch=legend.pch,col=colo,legend=lev,cex=1)
      }
    invisible(out)
  }
