\name{rotonto}
\alias{rotonto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{rotates and translates one matrix onto an other
}
\description{
rotate a matrix onto an other without loosing information about the location of the targetmatrix
}
\usage{
rotonto(x, y, scaling = FALSE, signref = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{k x m matrix to be rotated onto (targetmatrix)
%%     ~~Describe \code{x} here~~
}
  \item{y}{k x m matrix which will be rotated (reference matrix)
%%     ~~Describe \code{y} here~~
}
  \item{scaling}{logical: scale matrix to minimize sums of squares
%%     ~~Describe \code{scaling} here~~
}
  \item{signref}{logical: report if reflections were involved in the rotation
%%     ~~Describe \code{signref} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
  If it is a LIST, use
  	\item{yrot }{rotated and translated matrix}
 	 \item{Y }{centred reference matrix }
	\item{X }{centred target matrix }
	\item{trans }{vector between original position of target and centered reference (during rotation process) }
	\item{transy }{vector between original position of reference and centered reference (during rotation process) }	
	\item{gamm }{rotation matrix }	
	
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

function (x, y, scaling = FALSE, signref = T) 
{
    m <- dim(x)[2]
    X <- apply(x, 2, scale, scale = F)
    Y <- apply(y, 2, scale, scale = F)
    XY <- crossprod(X, Y)
    sv1 <- svd(XY)
    gamm <- (sv1$v) \%*\% t(sv1$u)
    if (sign(det(gamm)) < 1 && signref == T) {
        cat("reflection involved")
    }
    trans <- x[1, ] - X[1, ]
    transy <- y[1, ] - Y[1, ]
    sig <- sign(det(XY))
    del <- sv1$d
    ctrace <- function(MAT) sum(diag(crossprod(MAT)))
    bet <- sum(del)/ctrace(Y)
    if (scaling == TRUE) {
        yrot <- bet * Y \%*\% gamm
    }
    else {
        yrot <- Y \%*\% gamm
    }
    yrot <- t(apply(yrot, 1, function(x) {
        x + trans
    }))
    return(list(yrot = yrot, Y = Y, X = X, trans = trans, transy = transy, 
        gamm = gamm))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
