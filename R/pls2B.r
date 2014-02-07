#' Two-Block partial least square regression.
#' 
#' Performs a Two-Block PLS on two sets of data and assesses the significance
#' of each score by permutation testing
#' 
#' The Two-Block PLS tries to find those linear combinations in each block
#' maximising the covariance between blocks. The significance of each linear
#' combination is assessed by comparing the singular value to those obtained
#' from permuted blocks. If both blocks contain landmarks superimposed
#' TOGETHER, the option \code{same.config=TRUE} requests superimposition of the
#' permuted configurations (i.e. where the the landmarks of block \code{x} are
#' replaced by corresponding landmarks of other specimen.
#' 
#' @param y array containing superimposed landmark data of the first block.
#' Matrices are also allowed but the option 'same.config' will not work.
#' @param x array containing superimposed landmark data second block.Matrices
#' are also allowed but the option 'same.config' will not work.
#' @param tol threshold for discarding singular values.
#' @param same.config logical: if \code{TRUE} each permutation includes new
#' superimposition of permuted landmarks. This is necessary if both blocks
#' originate from landmarks that are superimposed together.
#' @param rounds rounds of permutation testing.
#' @param mc.cores integer: determines how many cores to use for the
#' computation. The default is autodetect. But in case, it doesn't work as
#' expected cores can be set manually. Parallel processing is disabled on
#' Windows due to occasional errors.
#' @return
#' \item{svd }{singular value decomposition (see \code{\link{svd}}) of the
#' 'common' covariance block}
#' \item{Xscores }{PLS-scores of x}
#' \item{Yscores }{PLS-scores of y}
#' \item{CoVar }{Dataframe containing singular values, explained
#' covariation, correlation coeffictient between PLS-scores and p-values}
#' @author Stefan Schlager
#' @seealso \code{\link{svd}}
#' @references Rohlf FJ, Corti M. 2000. Use of two-block partial least-squares
#' to study covariation in shape. Systematic Biology 49:740-753.
#' @keywords dynamic
#' @examples
#' 
#' library(shapes)
#' ### very arbitrary test:
#' ### check if first 4 landmarks covaries with the second 4
#' proc <- procSym(gorf.dat)
#' ## we do only 50 rounds to minimize computation time
#' \dontrun{#same.config takes too long for CRAN check
#' pls1 <- pls2B(proc$rotated[1:4,,],proc$rotated[5:8,,],
#'               same.config=TRUE,rounds=50,mc.cores=2)
#' }
#' pls1 <- pls2B(proc$rotated[1:4,,],proc$rotated[5:8,,],
#'               same.config=FALSE,rounds=50,mc.cores=2)
#' pls1$CoVar
#' layout(matrix(1:4,2,2,byrow=TRUE))
#' for(i in 1:4)
#'  plot(pls1$Xscores[,i]~pls1$Yscores[,i])
#' 
#' 
#' 
#' 
#' @export
pls2B <- function(x, y, tol=1e-12, same.config=FALSE, rounds=0, mc.cores=parallel::detectCores())
  {
    landmarks <- FALSE
    xorig <- x
    yorig <- y
    win <- FALSE
    if(.Platform$OS.type == "windows")
      win <- TRUE
    else
      registerDoParallel(cores=mc.cores)### register parallel backend
    
    if (length(dim(x)) == 3) {
        landmarks <- TRUE
        x <- vecx(x)
    }
    if (length(dim(y)) == 3)
      y <- vecx(y)
    else
      landmarks <- FALSE
      
    xdim <- dim(x)
    ydim <- dim(y)

    if (same.config && !landmarks)
      warning("the option same.config requires landmark array as input")
    
    
    cova <- cov(cbind(x,y))
    svd.cova <- svd(cova[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])

    svs <- svd.cova$d
    svs <- svs/sum(svs)
    svs <- svs[which(svs > 0.001)]

    covas <- svs*100
    l.covas <- length(covas)
    z1 <- x%*%svd.cova$u[,1:l.covas] #pls scores of x
    z2 <- y%*%svd.cova$v[,1:l.covas] #pls scores of y
    
### calculate correlations between pls scores
    cors <- 0
    for(i in 1:length(covas))
        cors[i] <- cor(z1[,i],z2[,i])
     

### Permutation testing
    permupls <- function(i)
      {
        x.sample <- sample(1:xdim[1])
        y.sample <- sample(x.sample)
        if (same.config && landmarks) {
            tmparr <- .bindArr2(xorig[,,x.sample],yorig[,,y.sample],along=1)
            tmpproc <- ProcGPA(tmparr,silent=TRUE)
            x1 <- vecx(tmpproc$rotated[1:dim(xorig)[1],,])
            y1 <- vecx(tmpproc$rotated[1:dim(yorig)[1],,])
        } else {
            x1 <- x
            y1 <- y
        }
        cova.tmp <- cov(cbind(x1[x.sample,],y1[y.sample,]))
        svd.cova.tmp <- svd(cova.tmp[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])
        svs.tmp <- svd.cova.tmp$d
        return(svs.tmp[1:l.covas])
      }
    p.values <- rep(NA,l.covas)
    if (rounds > 0) {
        if (win)
            permuscores <- foreach(i = 1:rounds, .combine = cbind) %do% permupls(i)
        else
            permuscores <- foreach(i = 1:rounds, .combine = cbind) %dopar% permupls(i)
        
        p.val <- function(x,rand.x)
            {
            p.value <- length(which(rand.x >= x))
            
            if (p.value > 0)
                p.value <- p.value/rounds
            else
                p.value <- 1/rounds
            gc()
            return(p.value)
          }
        
        for (i in 1:l.covas)
            p.values[i] <- p.val(svd.cova$d[i],permuscores[i,])
          
      }
### create covariance table
    Cova <- data.frame(svd.cova$d[1:l.covas],covas,cors,p.values)
    colnames(Cova) <- c("singular value","% total covar.","Corr. coefficient", "p-value")
    out <- list(svd=svd.cova,Xscores=z1,Yscores=z2,CoVar=Cova)
    return(out)
  }
