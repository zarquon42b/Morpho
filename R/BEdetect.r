BEdetect <- function(x,y)
  {

    L <- CreateL(x)
    coeff <- L$Lsubk%*%y
    coeff <- rbind(coeff,matrix(0,4,3))
    dif.be <- fx(x,y,coeff)
    dif.be <- dif.be^2
    dif.be <- apply(dif.be,1,sum)
    return(dif.be)
  }
