##' Create Matrices necessary for Thin-Plate Spline
#' 
#' Create (Bending Engergy) Matrices necessary for Thin-Plate Spline, and
#' sliding of Semilandmarks
#' 
#' 
#' @param matrix k x 3 or k x 2 matrix containing landmark coordinates.
#' @param lambda numeric: regularization factor
#' @param output character vector: select which matrices to create. Can be a vector containing any combination of the strings: \code{"K", "L","Linv","Lsubk", "Lsubk3"}.
#' @param threads threads to be used for parallel execution calculating K.
#' sliding of semilandmarks.
#' @return depending on the choices in \code{output}:
#' \item{L }{Matrix K as specified in Bookstein (1989)}
#' \item{L }{Matrix L as specified in Bookstein (1989)}
#' \item{Linv }{Inverse of matrix L as specified in Bookstein (1989)}
#' \item{Lsubk }{uper left k x k submatrix of \code{Linv}}
#' \item{Lsubk3 }{Matrix used for sliding in \code{\link{slider3d}} and \code{\link{relaxLM}}}.
#' @note This function is not intended to be called directly - except for
#' playing around to grasp the mechansims of the Thin-Plate Spline.
#' @seealso \code{\link{tps3d}}
#' @references Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005.
#' Semilandmarks in Three Dimensions, in Modern Morphometrics in Physical
#' Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer
#' Academic/Plenum Publishers.
#' 
#' Bookstein FL. 1989. Principal Warps: Thin-plate splines and the
#' decomposition of deformations. IEEE Transactions on pattern analysis and
#' machine intelligence 11(6).
#' 
#' @examples
#' 
#' data(boneData)
#' L <- CreateL(boneLM[,,1])
#' ## calculate Bending energy between first and second specimen:
#' be <- t(boneLM[,,2])%*%L$Lsubk%*%boneLM[,,2]
#' ## calculate Frobenius norm 
#' sqrt(sum(be^2))
#' ## the amount is dependant on on the squared scaling factor
#' # scale landmarks by factor 5 and compute bending energy matrix
#' be2 <- t(boneLM[,,2]*5)%*%L$Lsubk%*%(boneLM[,,2]*5)
#' sqrt(sum(be2^2)) # exactly 25 times the result from above
#' ## also this value is not symmetric:
#' L2 <- CreateL(boneLM[,,2])
#' be3 <- t(boneLM[,,1])%*%L2$Lsubk%*%boneLM[,,1]
#' sqrt(sum(be3^2))
#' @importFrom Matrix bdiag
#' @export
CreateL <- function(matrix,lambda=1e-8, output=c("K","L","Linv","Lsubk", "Lsubk3"),threads=1) {
    if (ncol(matrix) %in%  c(2,3)) {
        out <- list()
        k <- nrow(matrix)
        m <- ncol(matrix)
        Q <- cbind(1,matrix)
                                        #O <- matrix(0,4,4)
        if (!is.matrix(matrix) || !is.numeric(matrix))
            stop("matrix must be a numeric matrix")
        K <- .Call("createL",matrix,threads)
        L <- matrix(0,k+m+1,k+m+1)
        if (lambda !=0 )
            diag(K) <- lambda
                                        #K <- forceSymmetric(K)
        if ("K" %in% output)
            out$K <- K
        if (length(grep("L",output))) {
            L[1:k,1:k] <- K
            L[(k+1):(k+m+1),1:k] <- t(Q)
            L[1:k,(k+1):(k+m+1)] <- Q
            L <- forceSymmetric(L)
        }
        if ("L" %in% output)
            out$L <- L
        if ("Linv" %in% output || "Lsubk" %in% output || "Lsubk3" %in% output) {
            L1 <- try(solve(L),silent=TRUE)
            if (class(L1)=="try-error") {
                cat("CreateL: singular matrix: general inverse will be used.\n")
                L1 <- armaGinv(as.matrix(L))		
            }
            if ("Linv" %in% output)
                out$Linv <- L1
            if ("Lsubk" %in% output || "Lsubk3" %in% output)
                Lsubk <- forceSymmetric(L1[1:k,1:k])
            if ("Lsubk" %in% output)
                out$Lsubk <- Lsubk
            
            if ("Lsubk3" %in% output) {
                if (m == 3)
                    Lsubk3 <- forceSymmetric(bdiag(Lsubk,Lsubk,Lsubk))
                else
                    Lsubk3 <- forceSymmetric(bdiag(Lsubk,Lsubk))
                out$Lsubk3 <- Lsubk3
            }
        }
        
        return(out)
    } else
        stop("only works for matrices with 2 or 3 columns")
}

