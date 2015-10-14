#' compute the area of an n-dimensional hypersphere
#'
#' compute the area of an n-dimensional hypersphere
#'
#' @param n dimensionality of space the hypersphere is embedded in (e.g.3 for a 3D-sphere)
#' @param r radius of the sphere
#' @return
#' returns the area
#' @examples
#' areaSphere(2) #gives us the circumference of a circle of radius 1
#' @export
areaSphere <- function(n,r=1) {
    nom <- 2*pi^(n/2)*r^(n-1)
    denom <- gamma(n/2)
    return(nom/denom)
}

#' compute the area of an n-dimensional hypersphere cap
#' 
#' compute the area of an n-dimensional hypersphere cap
#' @param n  dimensionality of space the hypersphere is embedded in (e.g.3 for a 3D-sphere)
#' @param phi angle between vectors defining the cone
#' @param r radius of the sphere
#' @return
#' returns the area of the hypersphere cap
#' @examples
#' areaSpherePart(2,pi/2) # covers half the area of a circle
#' @export
areaSpherePart <- function(n,phi,r=1) {
    A <- areaSphere(n,r)
    ib <- pbeta(sin(phi)^2,(n-1)/2,0.5)
    Apart <- 0.5*ib*A
    if (phi > pi/2) {
        Apart <- A-Apart
    }
    return(Apart)
}

#' Test whether the direction of two vectors is similar
#'
#' Test whether the direction of two vectors is similar
#' @param x vector
#' @param y vector
#' @return
#' a list with
#' \item{angle}{angle between vectors}
#' \item{p.value}{p-value for the probability that the angle between two random vectors is smaller or equal to the one calculatted from x and y}
#' @details
#' Under the assumption of all (normalized) n-vectors being represented by an n-dimensional hypersphere, the probability of the angle between two vectors is <= the measured values can be estimated as the area of a cap defined by that angle and divided by the hypersphere's complete surface area.
#' @references
#' S. Li , 2011. Concise Formulas for the Area and Volume of a Hyperspherical Cap. Asian Journal of Mathematics & Statistics, 4: 66-70.
#' @examples
#' x <- c(1,0); y <- c(1,1) # for a circle this should give us p = 0.25 as the angle between vectors
#' ##  is pi/4 and for any vector the segment +-pi/4 covers a quarter of the circle 
#' angleTest(x,y)
#' @export
angleTest <- function(x,y) {
    if (length(x) != length(y))
        stop("vectors must be of the same length")
    phi <- angle.calc(x,y)
    n <- length(x)
    pval <- anglePvalue(phi,n)
    return(list(angle=phi,p.value=pval))
}

# returns the pvalue 
anglePvalue <- function(phi,n) {
    A <- areaSphere(n)
    Apart <- areaSpherePart(n=n,phi)
    out <- Apart/A
    return(out)
}
