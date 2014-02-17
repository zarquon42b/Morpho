#' calculate variance of a distribution stemming from prediction models
#' 
#' calculates a quotient of the overall varriance within a predicted
#' distribution to that from the original one. This function calculates a naive extension of the univariate R^2-value by
#' dividing the variance in the predicted dat by the variance of the original
#' data. No additional adjustments are made!!
#'
#' @param model a model of classes "lm" or "mvr" (from the package "pls")
#' @param ncomp How many latent variables to use (only for mvr models)
#' @param val use cross-vaildated predictions (only for mvr models)
#' @param \dots currently unused additional arguments.
#' @return 
#' returns the quotient.
#' @note The result is only!! a rough estimate of the variance explained by a
#' multivariate model. And the result can be misleading - especially when there
#' are many predictor variables involved. If one is interested in the value
#' each factor/covariate explains, we recommend a 50-50 MANOVA perfomed by the
#' R-package "ffmanova", which reports this value factor-wise.
#' @author Stefan Schlager
#' @references Langsrud O, Juergensen K, Ofstad R, Naes T. 2007. Analyzing
#' Designed Experiments with Multiple Responses Journal of Applied Statistics
#' 34:1275-1296.
#' 
#' @examples
#' 
#' lm1 <- lm(as.matrix(iris[,1:4]) ~ iris[,5])
#' exVar(lm1)
#' @rdname exVar
#' @export
exVar <- function(model,...)UseMethod("exVar")

#' @rdname exVar
#' @method exVar lm
#' @export
exVar.lm <- function(model,...)
  {
    exVar <- sum(eigen(cov(predict(model)))$values)/sum(eigen(cov(predict(model)+model$residuals))$values)
    return(exVar)
  }

#' @rdname exVar
#' @method exVar mvr
#' @export
exVar.mvr <- function(model,ncomp,val=FALSE,...)
  {
    if (!val || is.null(model$validation))
        exVar <- sum(Re(eigen(cov(model$fitted.values[,,ncomp]))$values))/sum(Re(eigen(cov(model$fitted.values[,,ncomp]+model$residuals[,,ncomp]),symmetric = T)$values))
    else
      {
        exVar <- sum(Re(eigen(cov(model$validation$pred[,,ncomp]))$values))/sum(Re(eigen(cov(model$fitted.values[,,ncomp]+model$residuals[,,ncomp]),symmetric = T)$values))
    }
    return(exVar)
  }
