predictShape.lm <- function(fit, datamod, PC, mshape)
{
  if (!inherits(fit, "lm"))
      stop("provide linear model of class 'lm'!")
  dims <- dim(mshape)
  if (!missing(datamod)) {
      mat <- model.matrix(datamod)
      pred <- mat%*%fit$coefficients
      names <- as.matrix(model.frame(datamod))
      names <-  apply(names,1, paste,collapse= "_")
      names <- gsub(" ","",names)
  } else {
      pred <- predict(fit)
      names <- rownames(pred)
  }
  predPC <- t(PC%*%t(pred))
  
  if (dim(pred)[1] > 1) {
      out <- array(NA,dim=c(dims,dim(pred)[1]))
      for (i in 1:dim(out)[3])
          out[,,i] <- mshape+matrix(predPC[i,],dims[1],dims[2])
  } else {
      out <- mshape+matrix(predPC,dims[1],dims[2])
  }
  if (length(dim(out)) == 3)
      dimnames(out)[[3]] <- names
  rownames(pred) <- names
  return(list(predicted=out,predictedPC=pred))
}
