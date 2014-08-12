#' classify specimen based on between-group PCA or CVA
#'
#' classify specimen based on between-group PCA or CVA
#'
#' @param x result of groupPCA or CVA
#' @param cv logical: use cross-validated scores if available
#' @return
#' \item{class}{classification result}
#' \item{groups}{original grouping variable}
#'
#' for object of CVA, also the posterior probabilities are returned.
#' @rdname classify
#' @export
classify <- function(x,cv=TRUE) UseMethod("classify")

#' @rdname classify
#' @export
classify.bgPCA <- function(x,cv=TRUE) {

    if (!is.null(x$CV) && cv)
        CV <- x$CV
    else
        CV <- x$Scores
    classVec <- NULL
    GmeanCenter <- sweep(x$groupmeans,2,x$Grandmean)
    GmeanScores <- GmeanCenter%*%x$groupPCs
    for (i in 1:nrow(CV)) {
        tmpdist <- (sqrt(rowSums(sweep(GmeanScores,2,CV[i,])^2)))
        classVec[i] <- names(tmpdist)[which(tmpdist == min(tmpdist))]
    }
    classVec <- factor(classVec)
    out <- list(class=classVec,groups=x$groups)
    class(out) <- "classify"
    return(out)
}

#' @rdname classify
#' @export
classify.CVA <- function(x,cv=T) {

    if (!is.null(x$class) && cv) {
        out <- list(class=x$class,groups=x$groups,posterior=x$posterior)
        class(out) <- "classify"
        return(out)
    } else {
        if (length(dim(x$Grandm)) == 2) {
            x$Grandm <- as.vector(x$Grandm)
            x$groupmeans <- vecx(x$groupmeans)
        }
        CV <- x$CVscores
        classVec <- NULL
        GmeanCenter <- sweep(x$groupmeans,2,x$Grandm)
        GmeanScores <- GmeanCenter%*%x$CV
        classprobs <- NULL
        for (i in 1:nrow(CV)) {
            tmpdist <- (rowSums(sweep(GmeanScores,2,CV[i,])^2))
            post <- probpost(tmpdist,x$prior)
            classVec[i] <- names(tmpdist)[which(post == max(post))]
            classprobs <- rbind(classprobs,post)
        }
        names(classVec) <- rownames(classprobs) <- rownames(x$CVscores)
        colnames(classprobs) <- rownames(x$groupmeans)
        classVec <- factor(classVec)
        out <- list(class=classVec,groups=x$groups,posterior=classprobs)
        class(out) <- "classify"
        return(out)
    }
}
    
    
#' @export       
print.classify <- function(x,...) {
    cat(" classification result in frequencies\n")
    tab <- table(x$groups,x$class)
    probtab <- prop.table(tab,1)
    print(tab)
    cat("\n\n classification result in %\n")
    print(probtab)

    cat(paste0("\n\n overall classification accuracy: ",100*round(mean(diag(probtab)),digits=5)," %\n"))
}
    
