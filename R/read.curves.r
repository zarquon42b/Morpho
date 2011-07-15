read.curves <- function(x)
  {
    tmp <- as.numeric(readLines(x,n=1))
    print(tmp)
    tmp1 <-  readLines(x,n=tmp+3)
    tmp2 <- strsplit(tmp1,split=" ")
    chckentr <- unlist(lapply(tmp2,length))
    vertex <- which(chckentr==4)
    lv <- length(vertex)
    vertices <- matrix(as.numeric(unlist(tmp2[vertex])),lv,4,byrow=TRUE)

    n <- as.integer(vertices[,4])
    vertices <- vertices[,1:3]

    return(list(vertices=vertices,id=n))
  }
