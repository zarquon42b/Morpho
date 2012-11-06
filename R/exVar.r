exVar <- function(model,...)UseMethod("exVar")

exVar.lm <-function(model,...)
  {
    exVar <-  sum(eigen(cov(predict(model)))$values)/sum(eigen(cov(predict(model)+model$residuals))$values)
    return(exVar)
  }
exVar.mvr <-function(model,ncomp,val=FALSE,...)
  {
    if (!val || is.null(model$validation))
        exVar <-  sum(Re(eigen(cov(model$fitted.values[,,ncomp]))$values))/sum(Re(eigen(cov(model$fitted.values[,,ncomp]+model$residuals[,,ncomp]),symmetric = T)$values))
    else
      {
        exVar <-  sum(Re(eigen(cov(model$validation$pred[,,ncomp]))$values))/sum(Re(eigen(cov(model$fitted.values[,,ncomp]+model$residuals[,,ncomp]),symmetric = T)$values))
    }
    return(exVar)
  }
