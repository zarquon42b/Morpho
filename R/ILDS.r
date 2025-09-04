#' Compute Pairwise Landmark Distances
#'
#' Compute pairwise landmark distances from landmark array
#' @param A k x m x n array of landmarks
#' @return
#' returns a matrix with each row containing the pairwise distances between landmarks
#' @examples
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' ### compute ILDS from Procrustes aligned landmarks.
#' edma <- ILDS(proc$rotated)
#' @export
ILDS <- function(A) {
    
    n=dim(A)[3]; p=dim(A)[1]; k=dim(A)[2] # n; p; k
    Name=NA
    for(i in 1:(p-1)) {
        Name[(sum(((p-1):1)[1:i])-p+i+1):(sum(((p-1):1) [1:i]))]=paste(i,(i+1):p,sep="-")
    }
    ildn=(1:(p*(p-1)/2)) # ildn # number of ILD with p landmarks
    ildl=Name[ildn] # ildl # names of ILD

    E <- t(apply(A,3,function(x) x <- as.vector(dist(x))))
    colnames(E)=ildl
    return(E)
}


#' Compute R2 for Interlandmark Distances
#'
#' Compute R2 for Interlandmark Distances explaining between group differences
#' @param x array containing landmarks
#' @param groups vector containing group assignments or a numeric covariate. For groups with more than two levels, a pair of needs to be specified using \code{which}
#' @param R2tol numeric: upper percentile for ILD R2 in relation to factor or in case \code{autocluster=TRUE}, the minimum quantile allowed.
#' @param autocluster logical: if TRUE, the function is trying to find a cluster with the highest R2-values. In this case \code{R2tol} is used as the quantile the cluster is allowed to occupy.
#' @param gap logical: if TRUE, the largest gap in the R2 distribution is sought. If FALSE, a clustering procedure is used.
#' @param bg.rounds numeric: number of permutation rounds to assess between group differences
#' @param wg.rounds numeric: number of rounds to assess noise within groups by bootstrapping.
#' @param which integer (optional): in case the factor levels are > 2 this determins which factorlevels to use
#' @param reference matrix containing start config landmarks. If NULL, it will be computed as mean for group 1.
#' @param target matrix containing target config landmarks. If NULL, it will be computed as mean for group 2.
#' @param mc.cores integer: number of cores to use for permutation tests.
#' @param plot logical: if TRUE show graphical output of steps involved
#' @param silent logical: suppress console output
#' @param ... additional parameters for internal use only.

#' @importFrom Morpho vecx bindArr
#' @return
#' A list containing:
#' \item{relevantILDs}{containing landmark information with the highest R2-values}
#' \item{allR2}{vector with ILD specific R2-values, sorted decreasingly}
#' \item{reftarILDS}{matrix with columns containing ILDs for reference and target shapes}
#' \item{sampleILD}{matrix containing ILDs of entire sample}
#' \item{R2tol}{R2-threshold used}
#' \item{reference}{reference used}
#' \item{target}{target used}
#' \item{bg.test}{result from between-group testing}
#' \item{confR2}{confidence for relevant ILDs from bootstrapping}
#' \item{x}{Procrustes superimposed raw data (without scaling)}
#' 
#' @examples
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' ilds <- ILDSR2(proc$rotated,groups,plot=TRUE,wg.rounds=99,mc.cores=1)
#' if (interactive())
#' visualize(ilds)
#' ## use covariate
#' \dontrun{
#' ildsLM <- ILDSR2(proc$rotated,groups=proc$size,plot=TRUE,wg.rounds=99,mc.cores=1)
#' if (interactive())
#' visualize(ildsLM)
#' }
#'
#' ## 2D Case with size as predictor
#' require(shapes)
#' require(Morpho)
#' gor.dat <- bindArr(gorf.dat,gorm.dat,along=3)
#' procg <- procSym(gor.dat)
#' ildsg <- ILDSR2(procg$rotated,procg$size,plot=FALSE,bg.rounds=999,wg.rounds=99,
#'                 mc.cores=1,autocluster=TRUE,R2tol=.8)
#' @importFrom Morpho permudist arrMean3
#' @import graphics stats 
#' @export 
ILDSR2 <- function(x,groups,R2tol=.9,autocluster=TRUE,gap=TRUE,bg.rounds=999,wg.rounds=999,which=1:2,reference=NULL,target=NULL,mc.cores=1,plot=FALSE,silent=FALSE,...) {
    D <- dim(x)[2] ## get LM dimensionality
        
    
    bootstrap <- FALSE
    
    args <- list(...)
    if ("bootstrap" %in% names(args)) {
        bootstrap <- args$bootstrap
    }
    if (!bootstrap)
        x <- procSym(x,CSinit = F,scale=F)$rotated
    
    ## set up grouping variable
    regression <- FALSE
    if (!is.factor(groups) && is.numeric(groups)) {
        regression <- TRUE
        if (!silent)
            message("groups is interpreted as covariate. Using regression scheme.")
    }
    
    if (is.factor(groups)) {
        groups <- factor(groups)
        lev <- levels(groups)
    }
    
    ## factor case
    if (!regression) {
        old_groups <- groups
        if (!is.null(which)) {
            groups <- factor(groups[groups %in% lev[which]])
            lev <- levels(groups)
            x <- x[,,which(old_groups %in% lev)]
          
        }
        ng <- length(lev)

        if (ng < 2)
            stop("provide at least two groups")
        if (length(groups) != dim(x)[3])
            warning("group affinity and sample size not corresponding!")

        ## create reference and target
        if (is.null(reference))
            reference <- arrMean3(x[,,groups==lev[1]])
        
        if (is.null(target))
            target <- arrMean3(x[,,groups==lev[2]])

        twosh <- bindArr(reference,target,along=3)
        xorig <- x
        x <- vecx(x,byrow = T)
    } else { ## regression case
        xorig <- x
        x <- vecx(x,byrow = T)
        lmod <- lm(x~groups)
        estquant <- quantile(groups,probs=c(.1,.9))
        tarref <- predict(lmod,newdata=list(groups=estquant))
        twosh <- vecx(tarref,revert = T,lmdim = D,byrow = T)
        reference <- twosh[,,1]
        target <- twosh[,,2]
    }
    
    ild <- ILDS(xorig)
    allILD <- round(ild, digits=6)
    mindim <- ncol(allILD)
    
    E <- ILDS(twosh)
    twosh.ILD <- round(as.data.frame(t(E)), digits=6)
    colnames(twosh.ILD)=c("start","target")

    av.twosh.ILD <- apply(twosh.ILD,1,mean)

    
    ratios.twosh.ILD <- twosh.ILD$target/twosh.ILD$start
    names(ratios.twosh.ILD) <- rownames(twosh.ILD)
    reftarILDratios <- sort(ratios.twosh.ILD)
    reftarMeanILDratios <- av.twosh.ILD[names(reftarILDratios)]
    
    ## compute R2
   # if (!regression)
  #      all.R2 <- as.vector(cor(allILD, as.numeric(groups))^2)
  #  else {
        R2lm <- (lm(allILD~groups))
        tmplm <- summary(R2lm)
        all.R2 <- sapply(tmplm,function(x) x <- x$r.squared)
 #   }
    names(all.R2) <- colnames(allILD)

    if (autocluster && !bootstrap) {
        R2tolOrig <- R2tol
        if (!silent)
            message("Autoclustering enabled. Trying to find cluster with highest R2-values")
        if (gap)
            nR2 <- findGap(all.R2)
        else
            nR2 <- findClusters(all.R2)
        R2tol <- 1-nR2/length(all.R2)
        if (R2tol < R2tolOrig) {
            R2tol <- R2tolOrig
            warning(paste("Autoclustering failed. Using provided R2tol value:",R2tol))
            
        } else {
            if (!silent)
                message(paste0("R2tol set to: ",round(R2tol,digits=2)))
        }
        
    }
    
    all.R2sorted <- sort(all.R2, decreasing=TRUE) # R2 of ILDs compared to factor in total sample
    reftarMeanILD <- av.twosh.ILD[names(all.R2sorted)]

    ## combine all sample wide R2 stats in a named list
    ILDstats <- list(reftarMeanILD=reftarMeanILD,reftarILDratios=reftarILDratios,reftarMeanILDratios=reftarMeanILDratios)
    
    ILDsR2 <- round(all.R2sorted[which(all.R2sorted > stats::quantile(all.R2sorted, probs=R2tol))], digits=7)
    ILDsRatios <- round(ratios.twosh.ILD[names(ILDsR2)], digits=7) # finds the corresponding ILDs ratios
    ILDsRatiosAbsRank <- 1+length(reftarILDratios)-rank(sort(round(abs(1-reftarILDratios), digits=7)), ties.method="random")[names(ILDsR2)]
    outOf100.ILDsRatiosAbsRank <- round(ILDsRatiosAbsRank*100/ncol(allILD), digits=0)

    ## create output table
    o1 <- round(rbind(ILDsR2, ILDsRatios, ILDsRatiosAbsRank, outOf100.ILDsRatiosAbsRank),digits=2)
    out <- list(relevantILDs=o1,allR2=all.R2sorted,reftarILDS=twosh.ILD,sampleILD=allILD,R2tol=R2tol,reference=reference,target=target,bg.rounds=bg.rounds,wg.rounds=wg.rounds,ILDstats=ILDstats)
    ## extract names of relevant ILDS
    R2names <- colnames(o1)

    ## compute between group permutation testing
    if (bg.rounds > 0 && !regression) {
        bg.test <- permudist(x,groups,rounds=bg.rounds)
        if (!silent)
            colorPVal(bg.test$p.value,rounds=bg.rounds)
        
        out$bg.test <- bg.test
    }
    if (regression && !bootstrap) {
        pca <- Morpho::prcompfast(x)
        bad <- which(pca$sdev^2 < 1e-12)
        if (length(bad))
            pca <- pca$x[,-bad]
        else
            pca <- pca$x
        
        R2lmPCA <- (lm(pca~groups))
        pval <- anova(R2lmPCA)$"Pr(>F)"[2]
        if (!silent)
            colorPVal(pval,permu=FALSE)
        out$bg.test$p.value <- pval
    }

    ## bootstrapping
    if (wg.rounds > 0) {
        wg.boot <- parallel::mclapply(1:wg.rounds,function(x) x <- bootstrapILDSR2(xorig,groups,R2tol=R2tol,regression=regression),mc.cores = mc.cores)
        out$wg.boot <- wg.boot
        freqsR2 <- unlist(lapply(wg.boot,match,R2names))
        confR2 <- sapply(1:length(R2names),function(x) x <- length(which(freqsR2==x)))
        confR2 <- round((((confR2+1)/(wg.rounds+1))*100),digits=3)
        names(confR2) <- R2names
        out$confR2 <- confR2
        
        if (!silent)
            colorILDS(confR2,wg.rounds,R2=all.R2)
    }     
    out$x <- xorig
    class(out) <- "ILDSR2"

    if (plot) {
        plot(out)
    }
    
    return(out)
}

#' @export
print.ILDSR2 <- function(x,...) {
    print(x$relevantILDs)
    cat("\n")
    if (x$bg.rounds > 0) {
        colorPVal(x$bg.test$p.value,rounds=x$bg.rounds)
    } 

    if (x$wg.rounds > 0) {
        colorILDS(x$confR2,x$wg.rounds,R2=x$allR2)
    }
}


bootstrapILDSR2 <- function(x,groups,R2tol,regression=FALSE) {
     xtmp <- x
    if (!regression) {
       
        lev <- levels(groups)
        for (i in lev) {
            tmpgroup <- which(groups==i)
            xtmp[,,groups==i] <- x[,,sample(tmpgroup,size=length(tmpgroup),replace = TRUE)]
        }
    } else {
        mysample <- sample(length(groups),replace=T,size=dim(x)[3]*2)
        xtmp <- x[,,mysample]
        groups <- groups[mysample]
    }
     out <- colnames(ILDSR2(xtmp,groups=groups,bg.rounds=0,wg.rounds=0,plot=FALSE,R2tol=R2tol,silent=TRUE,bootstrap=TRUE)$relevantILDs)
    
    
}
colorPVal <- function(x,rounds=NULL,permu=TRUE) {
    if (x < 0.05)
        pvalcol <- crayon::green
    else
        pvalcol <- crayon::red
    if (permu)
        message(crayon::bold(paste0("P-value between groups (",rounds," rounds): ",pvalcol(format(x,scientific=T,digits=3)),"\n")))
    else
        message(crayon::bold(paste0("P-value of linear model shape ~ predictor: ",pvalcol(format(x,scientific=T,digits=3)),"\n")))

    
}

colorILDS <- function(x,rounds=NULL,R2=NULL) {
    cat(crayon::bold(paste0("Bootstrapped (",rounds," rounds) confidence of relevant ILDS:\n")))
    
    for (i in 1:length(x)) {
        cat(" ")
        itmp <- x[i]
        itmp.name <- names(itmp)
        itmpR2 <- NULL
        if (!is.null(R2))
            itmpR2 <- round(R2[itmp.name],digits=2)
        if (itmp > 75)
            colfun <- crayon::green
        else if (itmp > 50)
            colfun <-  crayon::yellow
        else
            colfun <-  crayon::red
        cat(colfun(paste0(crayon::bold("ILDS",names(itmp)),":\t",itmp,"%\tR2=",itmpR2)))
        cat("\n")
    }
}



#' Plot the ILDS with the relevant ILDS ighlighted
#'
#' Plot the ILDS with the relevant ILDS ighlighted
#' @param x output of function \code{\link{ILDSR2}}
#' @param ref logical: if TRUE, the reference shape defined in  \code{\link{ILDSR2}} will be plotted. Otherwise the target is used. If \code{ref=FALSE}, the contraction and exxpansion will also be inverted
#' @param relcol color of relevant ILDs
#' @param rescol color of "irrelevant" ILDs
#' @param lwd numeric: define line width. Relevant ILDs are displayed by \code{3*lwd}.
#' @param cex numeric: size of plot content
#' @param col define color of landmarks
#' @param pch define symbols used to plot landmarks in 2D plot.
#' @param contractcol vector of colors for shortening ILDs, associated with confidence. Must be of \code{length(contol)}.
#' @param expandcol vector of colors for expanding ILDs, associated with confidence. Must be of \code{length(contol)}.
#' @param conftol vector: set thresholds for confidence coloring
#' @param useconf logical: if TRUE, highlighting according to supported confidence of ILD is applied.
#' @param add logical: if TRUE, plot is added to an existin one.
#' @param plot.legend logical: if TRUE, a legend is added to the plots with information on the coloring scheme.
#' @param ngrid integer: if \code{ngrid > 0}, a TPS grid is shown to display the spatial deformation from reference to target shape.
#' @param gridcol color of TPS grid
#' @param gridlty line type for TPS grid
#' @param links integer vector or list containing multiple integer vectors or a k x 2 integer matrix where each row defines a single link: add information on how landmarks are linked (aka wireframe)
#' @param linkcol color of wireframe defined by \code{links}
#' @param magnify numeric: symmetrically increase differences between reference and target by that factor
#' @param lollipop logical: if TRUE, landmark displacement between reference and target is displayed as lollipop graph.
#' @param lollipopcol color for lollipop handles
#' 
#' @param ... additional parameters passed to  \code{\link{deformGrid2d}} /  \code{\link{deformGrid3d}}.
#' @examples
#' ## 3D Example
#' require(Morpho)
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' ilds <- ILDSR2(proc$rotated,groups,plot=FALSE,bg.rounds=0,wg.rounds=99)
#' visualize(ilds)
#'
#' ## 2D Example
#' require(shapes)
#' gor.dat <- bindArr(gorf.dat,gorm.dat,along=3)
#' sex <- factor(c(rep("f",30),rep("m",29)))
#' procg <- procSym(gor.dat)
#' ildsg <- ILDSR2(procg$rotated,sex,plot=FALSE,bg.rounds=0,wg.rounds=99,R2tol=.9)
#' visualize(ildsg,cex=2,pch=19)
#'
#' ## use custom color and thresholds
#' visualize(ildsg,cex=2,pch=19,confcol=rainbow(5),contractcol=c(0.9,0.6,0.4))
#' @importFrom Morpho deformGrid2d deformGrid3d
#' @importFrom rgl text3d
#' @rdname visualize
#' @export
visualize <- function(x,...) UseMethod("visualize")

#' @export
#' @rdname visualize
#' @method visualize ILDSR2
visualize.ILDSR2 <- function(x,ref=TRUE,relcol="red",rescol="black",lwd=1,cex=2,col="red",pch=19,contractcol=c("red","orange"),expandcol=c("blue","cyan"),conftol=c(75,50),useconf=TRUE,add=FALSE,plot.legend=FALSE,ngrid=0,gridcol="grey90",gridlty=3,links=NULL,linkcol="grey50",lollipop=TRUE,lollipopcol="black", magnify=1,...) {
    if (!inherits(x, "ILDSR2")) 
        stop("please provide object of class 'ILDSR2'")
    reftarILDS <- x$reftarILDS
    rn <- rownames(reftarILDS)
    if (!is.null(links)) {
        if (is.matrix(links))
            links <- lapply(1:nrow(links),function(x) x <- links[x,])
    }
    pairing <- (matrix(as.integer(unlist(strsplit(rn,split = "-"))),length(rn),2,byrow=T))
    if (ref) {
        reference <- x$reference
        target <- x$target
    } else {
        reference <- x$target
        target <- x$reference
        x$ILDstats$reftarILDratios <- 1/x$ILDstats$reftarILDratios
    }
    if (magnify < 0)
        x$ILDstats$reftarILDratios <- 1/x$ILDstats$reftarILDratios
    if (magnify != 1) {
        mymean <- (reference+target)/2
        mydiff <- reference-mymean
        reference <- magnify*mydiff+mymean
        target <- -magnify*mydiff+mymean
    }
    ref0 <- reference[pairing[,1],]
    ref1 <- reference[pairing[,2],]
    tar0 <- target[pairing[,1],]
    tar1 <- target[pairing[,2],]
    if (!lollipop) {
        tar0 <- ref0
        tar1 <- ref1
    }
    D3 <- FALSE
    if (ncol(reference)==3) {
        D3 <- TRUE
        mydeform <- deformGrid3d
    } else
        mydeform <- deformGrid2d

    mydeform(reference,target,lines=FALSE,lwd=0,show=1,cex2=0,cex1=cex,col1=col,pch=pch,add=add,...)
    if (ngrid > 0)
        mydeform(reference,target,lines=F,lwd=0,show=1,cex2=0,cex1=0,add=TRUE,ngrid=ngrid,gridcol=gridcol,gridlty=gridlty,,...)
    if (is.null(x$confR2) || ! useconf) {
        highlight <- colnames(x$relevantILDs)
        if (!is.null(highlight)) {
            hm <- match(highlight,rn)
            
            mydeform(tar0[hm,,drop=FALSE],tar1[hm,,drop=FALSE],add=T,lcol = relcol,lwd=lwd*3,show=1,cex2=0,cex1=0,lty=1,...)
        }
    } else {
        leg.txt <- paste("Contraction - Conf. > ",conftol)
        leg.txt <- c(leg.txt,paste("Expansion - Conf. > ",conftol))
        highlight <- names(x$confR2)
        hm <- match(highlight,rn)
        myinterval <- getInterval(x$confR2,conftol)
        for (i in 1:(length(conftol))) {
            tmp <- which(myinterval == i)
            if (length(tmp)) {
                hmtmp <- match(highlight[tmp],rn)
                expand <- which(x$ILDstats$reftarILDratios[highlight[tmp]] > 1)
                contract <- which(x$ILDstats$reftarILDratios[highlight[tmp]] <= 1)
                if (length(expand))
                    mydeform(tar0[hmtmp[expand],,drop=FALSE],tar1[hmtmp[expand],,drop=FALSE],add=T,lcol = expandcol[i] ,lwd=lwd*3,show=1,cex2=0,cex1=0,lty=1,...)
                if(length(contract))
                    mydeform(tar0[hmtmp[contract],,drop=FALSE],tar1[hmtmp[contract],,drop=FALSE],add=T,lcol = contractcol[i] ,lwd=lwd*3,show=1,cex2=0,cex1=0,lty=1,...)
                
            }
        }
    }
   

     if (!is.null(links)) {
            if (lollipop)
                lineplot(target,point=links,lwd=lwd,add=TRUE,col=linkcol)
            else
                lineplot(reference,point=links,lwd=lwd,add=TRUE,col=linkcol)
        }
    
    if (D3) {
        if (lollipop)
              mydeform(reference,target,lines=T,lwd=lwd*2,show=1,cex2=0,cex1=0,add=TRUE,ngrid=0,lcol=lollipopcol,...)
         rgl::texts3d(reference,texts = 1:nrow(reference),adj=1.5,...)
         
         if (!is.null(x$confR2) && useconf && plot.legend) {
            if (interactive())
               answer <- readline("Please resize 3D window before legend is plotted and press any key")
             rgl::legend3d("topleft",leg.txt,col=c(contractcol,expandcol),title = "Confidence",lty=1,lwd=3)
         }
    }
    else {
        
        if (!lollipop)
            mydeform(reference,reference,lines=F,lwd=0,show=1,cex2=0,cex1=cex,col1=col,pch=pch,add=T,...)
        else
            mydeform(reference,target,lines=T,lwd=lwd*2,show=1,cex2=0,cex1=cex,col1=col,pch=pch,add=T,lty=1,lcol=lollipopcol,...)

        text(reference,adj=1,offset=1,cex=cex,...)
        ## ranges <- diff(range(reference[,1]))
        if (!is.null(x$confR2) && useconf )
            legend("topleft",bg="transparent",leg.txt,col=c(contractcol,expandcol),title = "Confidence",lty=1,lwd=3)
    }
} 

#' @rdname visualize
#' @export
visualise <- visualize


#' @rdname visualize
#' @export
visualise.ILDSR2 <- visualize.ILDSR2



#' Plot graphical report for ILDSR2
#'
#' Plot graphical report for ILDSR2
#' @param x output of \code{\link{ILDSR2}}
#' @param ... additonal parameter currently not used.
#' @export
plot.ILDSR2 <- function(x,...) {
    par(mfrow=c(2,2))

    plot(x$ILDstats$reftarMeanILDratios, x$ILDstats$reftarILDratios, main="ILD Ratio Variabilty vs. ILDs Values", xlab="Average of Start & Target ILD", ylab="Target/Start ILD Ratio")
    abline(a=1, b=0, col="grey", lwd=3, lty=1) 

    plot(x$ILDstats$reftarMeanILD,x$allR2, main="R2-Values vs. ILD Values", xlab="Average of Start & Target ILD", ylab="R2 for Sample ILDs vs Predictor")
    abline(a=quantile(x$allR2, probs=x$R2tol), b=0, col="grey", lwd=3, lty=1)

    hist(x$ILDstats$reftarILDratios, breaks=sqrt(length(x$ILDstats$reftarILDratios)), prob=TRUE, main="Disribution of Target/Start ILD ratios",xlab="Target/Start ILD Ratios")
    lines(density(x$ILDstats$reftarILDratios), col="red")

    hist(x$allR2, breaks=sqrt(length(x$allR2)), prob=TRUE, main=" R2-Value Distribution",xlab="R2-Values")
    lines(density(x$allR2), col="red")
    par(mfrow=c(1,1))
}


getInterval <- function(x,intervals) {    
    tmp <- sapply(x,function(x) x > intervals)
    tmp <- rbind(tmp,TRUE)
    chk <- apply(tmp,2,function(x) x <- which(x),simplify = FALSE)
    chk <- sapply(chk,min)
    return(chk)
    
}

#' @import mclust
findClusters <- function(x) {
    hc0 <- mclust::hc(x,c("E","V"))
    myclust <- mclust::Mclust(x,G=1:10,verbose = FALSE,modelNames = c("E","V"),initialization = list(hcPairs=hc0))
    m.best <- dim(myclust$z)[2]

    hc <- hclust(dist(x),method="ward.D2")
    tree <- cutree(hc,k=m.best)
    clustermeans <- sapply(min(tree):max(tree),function(y) y <- mean(x[tree=y]))
    bestC <- which.max(clustermeans)
    nbest <- length(which(tree == bestC))
    return(nbest)
                    
    
    
}

findGap <- function(x) {
    nl <- length(x)
    xsort <- sort(x,decreasing=T)
    xdiff <- xsort[1:(nl-1)]-xsort[2:nl]
    lgap <- which.max(xdiff)
    nbest <- length(which(xsort >= xsort[lgap]))
    return(nbest)
    
}
