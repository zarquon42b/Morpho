\name{read.lmdta}
\alias{read.lmdta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{read dta files
}
\description{
reads .dta files created by the software Landmark http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\usage{
read.lmdta(file = "x")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{a dta file 
%%     ~~Describe \code{file} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{arr}{array containing landmarks dimnames will be Information of ID and landmark names specified in Landmark}
 \item{info }{Information extracted from the header of the dta file}
 \item{idnames }{character vector containing the names of the individuals as specified in the dta file}
%% ...
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
