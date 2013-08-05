pred.PCmult <- function(fit,datamod,PC,mshape)
{
  
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat%*%fit$coefficients
  predPC <- t(PC%*%t(pred))
  if(dim(mat)[1] > 1)
    {out <- array(NA,dim=c(dims,dim(mat)[1]))
     for (i in 1:dim(out)[3])
       {out[,,i] <- mshape+matrix(predPC[i,],dims[1],dims[2])
                  
                  
      }
   }
     else
       {out <- mshape+matrix(predPC,dims[1],dims[2])
      }
	return(list(out=out,pred=pred))
}
pred.PC <- function (coeff, mod, PC, mshape) 
{
    dims <- dim(mshape)
    pred1 <- t(coeff) %*% mod
    predPC <- PC %*% pred1
    modell <- mshape + matrix(predPC, dims[1], dims[2])
    return(modell)
}
