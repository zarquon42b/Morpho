\name{read.csv.folder}
\alias{read.csv.folder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
batch import data from files
}
\description{
imports all data files contained in a specified folder.
}
\usage{
read.csv.folder(folder, x, y = 2:4, rownames = NULL, header = TRUE, dec = ".", sep = ";", pattern = "csv",addSpec=NULL, back = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder}{character: path to folder
%%     ~~Describe \code{folder} here~~
}
  \item{x}{either a vector specifiing which rows are to be imported, or character vector containing variable names to be sought for.
%%     ~~Describe \code{x} here~~
}
  \item{y}{a vector specifiing, which columns of the speradsheet ist to be imported.
%%     ~~Describe \code{y} here~~
}
  \item{rownames}{integer: specifies columns, where variable names are stored.
%%     ~~Describe \code{rownames} here~~
}
  \item{header}{logical : if spreadsheet contains header-row.
%%     ~~Describe \code{header} here~~
}
  \item{dec}{character: defines decimal sepearator.
%%     ~~Describe \code{dec} here~~
}
  \item{sep}{character: defines column seperator.
%%     ~~Describe \code{sep} here~~
}
  \item{pattern}{character: specify file format (e.g. csv).
%%     ~~Describe \code{pattern} here~~
}
  \item{addSpec}{character: add a custom specifier to the dimnames of the array.
%%     ~~Describe \code{pattern} here~~
}
  \item{back}{logical: where to place the specifier.
%%     ~~Describe \code{pattern} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{arr }{array containing imported data}
 \item{NAs }{vector containing position of observations with NAs}
 \item{NAs }{list: containing vectors containing information which LMs are missing in which observation}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Stefan Schlager
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

