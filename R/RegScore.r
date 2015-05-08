#' calulate regression scores for linear model
#'
#' calulate regression scores for linear model as specified in Drake & Klingenberg(2008)
#'
#' @param model linear model
#' @param x optional: matrix containing fitted data to be projected onto the regression lines. If omitted the model's fitted values will be used.
#' @return returns a n x m matrix containing the regression scores for each specimen.
#' @details the data are orthogonally projected onto the regression lines associated with each factor.
#' @section Warning: if \code{model} contains factors with more than 2 levels, R calculates one regression line per 2 factors. Check the \code{colnames} of the returned matrix to select the appropriate one. See examples for details.
#' @references Drake, AG. & Klingenberg, CP. The pace of morphological change: historical transformation of skull shape in St Bernard dogs. Proceedings of the Royal Society B: Biological Sciences, The Royal Society, 2008, 275, 71-76.
#' @examples
#' model <- lm(as.matrix(iris[,1:3]) ~ iris[,4])
#' rs <- RegScore(model)
#' plot(rs,iris[,4])
#'
#' ##now use a random subsample for model fitting
#' rand <- sample(nrow(iris))
#' x <- iris[rand[1:100],4]
#' newmod <- lm(as.matrix(iris[rand[1:100],1:3]) ~ x)
#' ##predict the rest of data and get their regression scores
#' rsPred <- RegScore(newmod,as.matrix(iris[rand[101:150],1:3]))
#' plot(rsPred,iris[rand[101:150],4])
#' \dontrun{
#' data(boneData)
#' proc <- procSym(boneLM)
#' pop.sex <- name2factor(boneLM,which=3:4) # generate a factor with 4 levels
#' lm.ps.size <- lm(proc$PCscores ~ pop.sex+proc$size)
#' rs <- RegScore(lm.ps.size)
#' colnames(rs) # in this case, the last column contains the regression
#' # scores associated with proc$size
#' ## validate by using a subsample for fitting
#' rand <- sample(dim(boneLM)[3])
#' lm.ps.size0 <- lm(proc$PCscores[rand[1:50],] ~ proc$size[rand[1:50]])
#' rs0 <- RegScore(lm.ps.size0,proc$PCscores[rand[-c(1:50)],] )
#' plot(rs0,proc$size[rand[-c(1:50)]])
#' }
#' 
#' @export

RegScore <- function(model,x=NULL) {
    if (is.null(x))
        x <- model$fitted.values+model$residuals
    coef <- coef(model)[-1,]
    if (!is.vector(coef))
        coef <- t(coef)

    if (is.matrix(coef))
        coefscale <- apply(coef,2,base::norm,"2")
    else
        coefscale <- base::norm(coef,"2")

    
    coef <- scale(coef,center=F,scale=coefscale)
    RegScores <- x%*%coef
    colnames(RegScores) <- colnames(coef)
    return(RegScores)
}
        
