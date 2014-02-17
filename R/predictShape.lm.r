#' Predict shapes based on linear models calculated from PCscores
#' 
#' Predict shapes based on linear models calculated from PCscores.
#' 
#' This function predicts the landmarks based on models calculated from
#' PCscores.
#' 
#' @param fit model of class \code{lm} where the PCscores are fitted onto
#' @param datamod a one-sided "model" formula, of the form \code{~ x1 + x2 +
#' ... + xk}, corresponding to the right hand term in the model used in
#' \code{fit}. If omitted, the predicted shapes of all specimen are calculated
#' based on the fitted values.
#' @param PC Matrix/vector containing Principal components (rotation matrix)
#' corresponding to PC-scores used in \code{fit}.
#' @param mshape matrix of the meanshape's landmarks by which the data was
#' centered before rotation in covariance eigenspace.
#' @return
#' \item{predicted }{array or matrix containing predicted landmark
#' coordinates}
#' \item{predictedPC }{matrix containing predicted PC-scores}
#' @section Warning: Make sure that the levels of the variables used in
#' \code{datamod} correspond exactly to those used in \code{fit}. Otherwise
#' model matrix will be calculated erroneous.
#' @seealso \code{\link{model.matrix}, \link{lm}, \link{formula}}
#' 
#' @examples
#' 
#' data(boneData)
#' proc <- procSym(boneLM)
#' pop <- name2factor(boneLM,which=3)
#' ##easy model with only one factor based on the first four PCs
#' fit <- lm(proc$PCscores[,1:4] ~ pop)
#' ## get shape for Europeans only
#' datamod <- ~as.factor(levels(pop))[2]
#' Eu <- predictShape.lm(fit,datamod, proc$PCs[,1:4],proc$mshape)
#' 
#' ## get shape for Europeans and Chinese
#' datamod <- ~as.factor(levels(pop))
#' pred <- predictShape.lm(fit,datamod, proc$PCs[,1:4],proc$mshape)
#' \dontrun{
#' deformGrid3d(pred$predicted[,,1], pred$predicted[,,2], ngrid = 0)
#' }
#' 
#' ## more complicated model
#' 
#' sex <- name2factor(boneLM,which=4)
#' fit <- lm(proc$PCscores[,1:4] ~ pop*sex)
#' ## predict female for chinese and European
#' datamod <- ~(as.factor(levels(pop))*rep(as.factor(levels(sex))[1],2))
#' pred <- predictShape.lm(fit,datamod, proc$PCs[,1:4],proc$mshape)
#' 
#' ## predict female and malefor chinese and European
#' popmod <- factor(c(rep("eu",2),rep("ch",2)))
#' sexmod <- rep(as.factor(levels(sex)),2)
#' datamod <- ~(popmod*sexmod)
#' pred <- predictShape.lm(fit,datamod, proc$PCs[,1:4],proc$mshape)
#' 
#' ## add some (randomly generated) numeric covariate
#' somevalue <- rnorm(80,sd=10)
#' fit <- lm(proc$PCscores[,1:4] ~ pop+somevalue)
#' probs <- quantile(somevalue, probs=c(0.05, 0.95))
#' ## make model for European at 5% and 95% quantile
#' popmod <- rep(factor(levels(pop))[2],2)
#' datamod <- ~(popmod+probs)
#' pred <- predictShape.lm(fit,datamod, proc$PCs[,1:4],proc$mshape)
#' 
#' 
#' @export
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
