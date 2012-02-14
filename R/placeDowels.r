placeDowels <- function(lm,mesh,ldowel,smooth=TRUE,render=TRUE,col=1,radius=1,meshcol=3,output=TRUE,fine=50)
  {

    ldo <- FALSE
    colvec <- FALSE
    projLM <- projRead(lm,mesh,smooth=smooth)
    if (length(ldowel) > 1)
      {
        ldo <- TRUE
      }
    
     if (length(col)>1)
       {
         colvec <- TRUE
       }
    
    dowels <- list()
### create dowels and render them if required
    for (i in 1:dim(lm)[1])
      {
        if (ldo)
          {
            ltmp <-ldowel[i]
          }
        else
          {
            ltmp=ldowel
          }
        dowels[[i]] <- cylinder(projLM$vb[1:3,i],projLM$normals[1:3,i],length=ltmp,radius=radius,fine=fine,adNormals=FALSE)
        if (render)
          {
           if (colvec)
             {
               coltmp <- col[i]
             }
           else
             {
               coltmp <- col
             }
            shade3d(dowels[[i]],col=coltmp)
         }

        
      }
    if (render) ## render mesh ##
      {
        shade3d(mesh,col=meshcol)
      }
    if (output)
      {
        return(dowels)
      }
  }
