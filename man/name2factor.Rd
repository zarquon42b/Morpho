\name{name2factor}
\alias{name2factor}
\alias{name2num}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{extract data from array names
}
\description{
extract data from array names and convert to factors or numbers}
\usage{
name2factor(x, sep = "_", which, collapse = sep, dif = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data, can be a three-dimensional array, a  matrix, a named list or a vector containing
    names to split
}
  \item{sep}{character by which to split the strings
%%     ~~Describe \code{sep} here~~
}
  \item{which}{integer or vector of integers, if more entries are
    selected, they will be concatenated by the string specified with the
    option 'collapse'.
%%     ~~Describe \code{which} here~~
}
  \item{collapse}{character by which to collapse data if two strings are
    to be concatenated
%%     ~~Describe \code{collapse} here~~
}
  \item{dif}{calculate difference if two fields containing numbers are selected.
}
}
\details{If an array is used as input, the data info is expected to be
  in the 3rd dimension, for a matrix, rownames are used.
%%  ~~ If necessary, more details than the description above ~~
}
\value{returns a vector containing factors or numbers
}

\author{
Stefan Schlager}



\examples{

data <- matrix(rnorm(200),100,2)
id <- paste("id",1:100,sep="")
pop <- c(rep("pop1",50),rep("pop2",50))
sex <- c(rep("male",50),rep("female",50))
age <- floor(rnorm(100,mean=50,sd=10))
rownames(data) <- paste(id,pop,sex,age,sep="_")
infos <- data.frame(pop=name2factor(data,which=2))
infos$age <- name2num(data,which=4)
infos$pop.sex <- name2factor(data,which=2:3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
