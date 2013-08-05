 angle.calc <- function(x,y)
 {  x <- as.vector(x)/sqrt(sum(x^2))
    y <- as.vector(y)/sqrt(sum(y^2))
    rho <- acos((sum((x-y)^2)-2)/-2)
    return(rho)
 }
