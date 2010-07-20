\name{render.obj}
\Rdversion{1.1}
\alias{render.obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rendering a 3D Mesh

}
\description{
renders an object in 3D, imported with read.obj and gives out an shapelist object.
}
\usage{
render.obj(obj, wire = "s", col = 1,plot = T, shapeout = T)

}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ object imported via "read.obj" command
%%     ~~Describe \code{obj} here~~
}
  \item{wire}{Option how the input is to be rendered if a graphical output is demanded. Otions are: "s"= surface, "w"=wireframe or "p"=cloud of points. 
%%     ~~Describe \code{wire} here~~
}
  \item{col}{color of the rendered visualization.
%%     ~~Describe \code{col} here~~
}
  
  \item{plot}{Logical: requests a 3D-output.
%%     ~~Describe \code{plot} here~~
}
  \item{shapeout}{Logical: if True, a shapelist3d object is generated.
%%     ~~Describe \code{shapeout} here~~
}
}
\details{ trimesh.obj is very much faster than render.obj but only works on triangular meshes
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{shape }{shapelist3d object}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{shapelist3d}},\code{\link{warp_obj}} ,\code{\link{trimesh.obj}},\code{\link{read.obj}} ,\code{\link{write.obj}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj, wire = "s", col = 1, size = size, plot = T, shapeout = T) 
{
    library(rgl)
    vert <- as.matrix(subset(obj, obj[, 1] == "v")[, 2:4])
    tri <- as.matrix(subset(obj, obj[, 1] == "f")[, 2:4])
    if (wire == "p") {
        points3d(vert, col = col, size = size)
    }
    else {
        meshlist <- list(0)
        for (i in 1:dim(tri)[1]) {
            meshlist[[i]] <- tmesh3d(c(vert[tri[i, 1], ], vert[tri[i, 
                2], ], vert[tri[i, 3], ]), indices = c(1:3), 
                homogeneous = F)
        }
        shape <- shapelist3d(meshlist, plot = F)
        if (wire == "w" && plot == TRUE) {
            wire3d(shape, color = col)
        }
        if (wire == "s" && plot == TRUE) {
            shade3d(shape, color = col)
        }
    }
    if (shapeout == T) {
        return(shape = shape)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

