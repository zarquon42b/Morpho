\name{CVA}
\Rdversion{1.1}
\alias{CVA}
\title{Canonical Variate Analysis
}
\description{
performs a Canonical Variate Analysis.
}
\usage{
CVA(dataarray, groups, weighting = TRUE, tolinv = 1e-10, plot = TRUE,
    rounds=0, cv=FALSE, mc.cores=detectCores())
}
\arguments{
  \item{dataarray}{Either a k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size. Or alternatively a n x m Matrix where n is the numeber of observations and m the number of variables (this can be PC scores for example)
}
  \item{groups}{a character/factor vector containgin grouping variable.
}
  \item{weighting}{Logical: Determines whether the between group covariance matrix is to be weighted according to group size.
}
  \item{tolinv}{Threshold for the eigenvalues of the pooled within-group-covariance matrix to be taken as zero - for calculating the general inverse of the pooled withing groups covariance matrix.
}
  \item{plot}{Logical: determins whether in the two-sample case a histogramm ist to be plotted.
}
 \item{rounds}{integer: number of permutations if a permutation test of
   the mahalanobis and Procrustes distance between group means is requested.If rounds = 0, no test is performed.
}
 \item{cv}{logical: requests a Jackknife Crossvalidation.
}
\item{mc.cores}{integer: how many cores of the Computer are allowed to be
    used. Default is use autodetection by using detectCores()
    from the parallel package.
}
}

\value{
        
\item{CV }{A matrix containing the Canonical Variates}
\item{CVscores }{A matrix containing the individual Canonical Variate scores}
\item{Grandm}{a vector or a matrix containing the Grand Mean (depending if the input is an array or a matrix)}
\item{groupmeans }{a matrix or an array containing the group means (depending if the input is an array or a matrix)}
\item{Var }{Variance explained by the Canonical Variates}
\item{CVvis }{Canonical Variates projected back into the original space - to be used for visualization purposes, for details see example below.}
\item{Dist }{Mahalanobis Distances between group means - if requested tested by permutation test if the input is an array it is assumed to be superimposed Landmark Data and Procrustes Distance will be calculated.}
\item{CVcv }{A matrix containing crossvalidated CV scores}
\item{mc.cores}{integer: determines how many cores to use for the
  computation. The default is to autodetect. But in case, it doesn't work as expected cores can be
  set manually.Parallel processing is disabled on
    Windows due to occasional errors.
}
}
\references{Cambell, N. A. & Atchley, W. R.. 1981 The Geometry of
  Canonical Variate Analysis: Syst. Zool., 30(3), 268-280.
  
Klingenberg, C. P. & Monteiro, L. R. 2005 Distances and directions in multidimensional shape spaces: implications for morphometric applications. Systematic Biology 54, 678-688.
}
\author{
Stefan Schlager
}


\seealso{
\code{\link{groupPCA}}
}
\examples{
## all examples are kindly provided by Marta Rufino

library(shapes)
# perform procrustes fit on raw data
alldat<-procSym(abind(gorf.dat,gorm.dat))
# create factors
groups<-as.factor(c(rep("female",30),rep("male",29)))
# perform CVA and test Mahalanobis distance
# between groups with permutation test by 100 rounds)            
cvall<-CVA(alldat$orpdata,groups,rounds=100,mc.cores=2)     

### Morpho CVA
data(iris)
vari=iris[,1:4]
facto=iris[,5]

#note that the function takes time, to estimate permutations.
cva.1=CVA(vari, groups=facto,mc.cores=2) 
# plot the CVA
plot(cva.1$CVscores, col=facto, pch=as.numeric(facto), typ="n",asp=1,
   xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"\%")),
   ylab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"\%")))
  
  text(cva.1$CVscores, as.character(facto), col=as.numeric(facto), cex=.7)

  # add chull (merge groups)
  for(jj in 1:length(levels(facto))){
        ii=levels(facto)[jj]
    kk=chull(cva.1$CVscores[facto==ii,1:2])
    lines(cva.1$CVscores[facto==ii,1][c(kk, kk[1])],
    cva.1$CVscores[facto==ii,2][c(kk, kk[1])], col=jj)
    }

  # add 80% ellipses
  require(car)
  for(ii in 1:length(levels(facto))){
    dataEllipse(cva.1$CVscores[facto==levels(facto)[ii],1],
    cva.1$CVscores[facto==levels(facto)[ii],2], 
                    add=TRUE,levels=.80, col=c(1:7)[ii])}

  # histogram per group
  require(lattice)
  histogram(~cva.1$CVscores[,1]|facto,
  layout=c(1,length(levels(facto))),
          xlab=paste("1st canonical axis", paste(round(cva.1$Var[1,2],1),"\%")))
  histogram(~cva.1$CVscores[,2]|facto, layout=c(1,length(levels(facto))),
          xlab=paste("2nd canonical axis", paste(round(cva.1$Var[2,2],1),"\%")))

  # plot Mahalahobis
  dendroS=hclust(cva.1$Dist$GroupdistMaha)
  dendroS$labels=levels(facto)
  par(mar=c(4,4.5,1,1))
  dendroS=as.dendrogram(dendroS)
  plot(dendroS, main='',sub='', xlab="Geographic areas",
          ylab='Mahalahobis distance')

 
   # Variance explained by the canonical roots:
   cva.1$Var
   # or plot it:
   barplot(cva.1$Var[,2])
 
 }
