#' classify specimen based on between-group PCA or CVA or typprobClass
#'
#' classify specimen based on between-group PCA, CVA or typprobClass
#'
#' @param x result of groupPCA, CVA or typprobClass
#' @param cv logical: use cross-validated scores if available
#' @return
#' \item{class}{classification result}
#' \item{groups}{original grouping variable}
#'
#' for object of CVA and typprob, also the posterior probabilities are returned.
#' @rdname classify
#' @seealso \code{\link{CVA},\link{groupPCA},  \link{typprobClass}}
#' @export
classify <- function(x,cv=TRUE) UseMethod("classify")

#' @rdname classify
#' @export
classify.bgPCA <- function(x,cv=TRUE) {

    if (!is.null(x$CV) && cv)
        CV <- x$CV
    else
        CV <- x$Scores
    if (is.null(x$CV)) 
        cv  <- FALSE
    classVec <- x$groups
    GmeanCenter <- sweep(x$groupmeans,2,x$Grandmean)
    GmeanScores <- GmeanCenter%*%x$groupPCs
    for (i in 1:nrow(CV)) {
        tmpdist <- (sqrt(rowSums(sweep(GmeanScores,2,CV[i,])^2)))
        classVec[i] <- names(tmpdist)[which(tmpdist == min(tmpdist))]
    }
    out <- list(class=classVec,groups=x$groups,cv=cv)
    class(out) <- "classify"
    return(out)
}

#' @rdname classify
#' @export
classify.CVA <- function(x,cv=T) {

    if (!is.null(x$class) && cv) {
        if (is.null(x$CVcv)) 
            cv  <- FALSE
        out <- list(class=x$class,groups=x$groups,posterior=x$posterior,cv=cv)
        class(out) <- "classify"
        return(out)
    } else {
        if (length(dim(x$Grandm)) == 2) {
            x$Grandm <- as.vector(x$Grandm)
            x$groupmeans <- vecx(x$groupmeans)
        }
        CV <- x$CVscores
        classVec <- x$groups
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
                                        #classVec <- factor(classVec)
        out <- list(class=classVec,groups=x$groups,posterior=classprobs,cv=cv)
        class(out) <- "classify"
        return(out)
    }
}

#' @rdname classify
#' @export
classify.typprob <- function(x,cv=TRUE) {
    out <- list()
    if (!x$cv)
        cv <- FALSE
    if (cv) {
        out$class <- x$groupaffinCV
        out$groups <- x$groups
        out$posterior <- x$probsCV
    } else {        
        out$class <- x$groupaffin
        out$groups <- x$groups
        out$posterior <- x$classprobs
    }
    out$cv=cv
    class(out) <- "classify"
    attributes(out) <- append(attributes(out),list(self=x$self))
    return(out)
}   

#' @export       
print.classify <- function(x,...) {
    if (is.null(attributes(x)$self))
        self <- TRUE
    else
        self <- attributes(x)$self

    out <- list()
    if (self) {
        if (x$cv)
            cat(" cross-validated classification results in frequencies\n")
        else 
            cat(" classification result in frequencies\n")
        tab <- table(x$groups,x$class)
        acc <- 0
        good <- which(as.character(x$groups) == as.character(x$class))
        acc <- 100*length(good)/length(x$groups)
        
        probtab <- prop.table(tab,1)*100
        print(tab)
        if (x$cv)
            cat("\n\n cross-validated classification result in %\n")
        else
            cat("\n\n classification result in %\n")
        print(probtab,digits=5)
        out$accuracy <- acc
        out$cv <- x$cv
        out$probtable <- probtab
        out$table <- tab
        cat(paste0("\n\n overall classification accuracy: ",round(acc,digits=5)," %\n"))
        invisible(out)
    } else {
        numbs <- tapply(x$class,x$class,length)
        cat("\n\n specimens have been classified as:\n\n ")
        print(numbs)
    }
                                        
}

