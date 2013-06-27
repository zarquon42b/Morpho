\name{ply2mesh}
\alias{ply2mesh}
\alias{obj2mesh}
\alias{file2mesh}
\title{
Import 3D surface mesh files
}
\description{
imports 3D mesh files and store them as an R .object of class mesh3d
}
\usage{
ply2mesh(filename, adnormals = TRUE,readnormals = FALSE,readcol=FALSE)
obj2mesh(filename, adnormals = TRUE)
file2mesh(filename, clean = TRUE,readcol=FALSE)
}
\arguments{
  \item{filename}{character: path to file
}
  \item{adnormals}{Logical: If the file does not contain normal information, they will be calculated in R: Can take some time.
}
\item{readnormals}{Logical: Import vertex normals (if available), although no face information is present.
}
\item{readcol}{Logical: Import vertex colors (if available).
}
\item{clean}{Logical: Delete dumpfiles.
}
}
\details{
file2mesh needs the command line tool ply2ascii from the Auxiliaries section in \url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/} installed.
}
\value{
\item{mesh }{list of class mesh3d - see rgl manual for further details, or a matrix containing vertex information or a list containing vertex and normal information}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
