#' calulate regression scores for linear model
#'
#' calulate regression scores for linear model as specified in Drake & Klingenberg(2008)
#'
#' @param model linear model
#' @param x optional: matrix containing new data to be projected onto the regression lines.
#' @return returns a n x m matrix containing the regression scores for each specimen.
#' @details the data is orthogonally projected onto the regression lines associated with each factor.
#' @references Drake, AG. & Klingenberg, CP. The pace of morphological change: historical transformation of skull shape in St Bernard dogs. Proceedings of the Royal Society B: Biological Sciences, The Royal Society, 2008, 275, 71-76.
#' @examples
#' model <- lm(as.matrix(iris[,1:2]) ~ iris[,3])
#' rs <- RegScore(model)
#' plot(rs,iris[,4])
#' @export RegScore

RegScore <- function(model,x=NULL) {
    if (is.null(x))
        x <- model$fitted.values+model$residuals
    coef <- coef(model)[-1,]
    if (!is.vector(coef))
        coef <- t(coef)
        RegScores <- x%*%coef
    colnames(RegScores) <- attr(model$terms,which="term.labels")
    return(RegScores)
}
        
