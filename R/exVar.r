exVar <-function(model)
  {
    exVar <-  sum(eigen(cov(predict(model)))$values)/sum(eigen(cov(predict(model)+model$residuals))$values)
    return(exVar)
  }
