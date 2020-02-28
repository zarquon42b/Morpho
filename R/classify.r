#' classify specimen based on between-group PCA or CVA or typprobClass
#'
#' classify specimen based on between-group PCA, CVA or typprobClass
#'
#' @param x result of groupPCA, CVA or typprobClass
#' @param cv logical: use cross-validated scores if available
#' @param newdata use new data to predict scores and evaluate group affinity
#' @param prior specify prior probability for CVA evaluation if NULL prior from CVA will be used. Be \code{m} your number of groups then to set the prior equally for all groups set \code{prior=rep(1,m)/m}.
#' @param ... currently not used
#' @return
#' \item{class}{classification result}
#' \item{groups}{original grouping variable, only available if \code{newdata=NULL}}
#' \item{posterior}{only for object of CVA and typprob, also the posterior probabilities are returned}
#' @rdname classify
#' @seealso \code{\link{CVA},\link{groupPCA},  \link{typprobClass}}
#' @export
classify <- function(x,cv=TRUE,...) UseMethod("classify")

#' @rdname classify
#' @export
classify.bgPCA <- function(x,cv=TRUE,newdata=NULL,...) {

    if (length(dim(x$groupmeans)) == 3) {
        x$groupmeans <- vecx(x$groupmeans)
        x$Grandmean <- c(x$Grandmean)
    }
    usenew <- FALSE
    if (is.null(newdata)) {
        if (!is.null(x$CV) && cv)
            CV <- x$CV
        else
            CV <- x$Scores
    } else {
        CV <- predict(x,newdata)
        usenew <- TRUE
    }
    if (is.null(x$CV) || usenew) 
        cv  <- FALSE

    classVec <- NULL
    GmeanCenter <- sweep(x$groupmeans,2,x$Grandmean)
    GmeanScores <- GmeanCenter%*%x$groupPCs
    for (i in 1:nrow(CV)) {
        tmpdist <- (sqrt(rowSums(sweep(GmeanScores,2,CV[i,])^2)))
        classVec[i] <- names(tmpdist)[which(tmpdist == min(tmpdist))]
    }
    
    outgroups <- x$groups
    
    if (usenew) {
        out <- list(class=classVec)
        attributes(out) <- append(attributes(out),list(self=FALSE))
    } else
        out <- list(class=classVec,groups=outgroups,cv=cv)
    class(out) <- "classify"
    return(out)
}

#' @rdname classify
#' @export
classify.CVA <- function(x,cv=T,newdata=NULL,prior=NULL,...) {
    
    if (!is.null(x$class) && cv && is.null(newdata)) {
        if (is.null(x$CVcv)) 
            cv  <- FALSE
        out <- list(class=x$class,groups=x$groups,posterior=x$posterior,cv=cv)
        class(out) <- "classify"
        return(out)
    } else {
        xorig <- x
        if (length(dim(x$Grandm)) == 2) {
            x$Grandm <- as.vector(x$Grandm)
            x$groupmeans <- vecx(x$groupmeans)
        }
        usenew <- FALSE
        if (is.null(newdata)) {
            CV <- x$CVscores
        }  else {
            CV <- predict(xorig,newdata=newdata)
            usenew <- TRUE
        }
        classVec <- NULL
        GmeanCenter <- sweep(x$groupmeans,2,x$Grandm)
        GmeanScores <- GmeanCenter%*%x$CV
        classprobs <- NULL
        if (is.null(prior))
            prior <- x$prior
        for (i in 1:nrow(CV)) {
            tmpdist <- (rowSums(sweep(GmeanScores,2,CV[i,])^2))
            post <- probpost(tmpdist,prior)
            classVec[i] <- names(tmpdist)[which(post == max(post))]
            classprobs <- rbind(classprobs,post)
        }
        if (!usenew)
            names(classVec) <- rownames(classprobs) <- rownames(x$CVscores)
        
        colnames(classprobs) <- rownames(x$groupmeans)
                                        #classVec <- factor(classVec)
        if (usenew) {
            out <- list(class=classVec,posterior=classprobs)
            attributes(out) <- append(attributes(out),list(self=FALSE))
        } else
            out <- list(class=classVec,groups=x$groups,posterior=classprobs,cv=cv)
        class(out) <- "classify"
        
        return(out)
    }
}

#' @rdname classify
#' @export
classify.typprob <- function(x,cv=TRUE,...) {
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
        ## Kappa statistics
        if (length(levels(x$groups)) == length(levels(x$class))) {
            n <- length(x$groups)
            p <- tapply(x$groups,x$groups,length)/n
            q <- tapply(x$class,x$class,length)/n
            expAccuracy <- sum(p*q)
            out$kappa = (out$accuracy/100 - expAccuracy) / (1 - expAccuracy)
            cat(paste0("\n Kappa statistic: ",round(out$kappa,digits=5),"\n"))
        }

        invisible(out)
    } else {
        numbs <- tapply(x$class,x$class,length)
        cat("\n\n specimens have been classified as:\n\n ")
        print(numbs)
    }
                                        
}

