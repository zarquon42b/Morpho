ray2mesh<-function(mesh1,tarmesh,tol=1,angmax=NULL,clean=TRUE,outname=NULL,readback=TRUE,inbound=FALSE,strict=FALSE,ignore.stdout=FALSE,mindist=FALSE)
{ 

  options <- NULL
  opt <- FALSE
  target <- "target.ply"
  if (inbound == TRUE)
    { opt <- TRUE
      options <- "--inbound" # set option to search along negative nromals first
    }
  if (strict == TRUE)
    {
      opt <- TRUE
      options <- paste(options,"--strict") #mark vertices that are not hit along rays
   }
  if (mindist == TRUE)
    {
      opt <- TRUE
      options <- paste(options,"--minray") #use closest point in and outward
    }
  if (!is.null(angmax))
    {
      opt <- TRUE
      options <- paste(options,"--angmax",angmax) #check for normals
    }
  if (opt)
    { options <- paste(" ",options,sep="")
    }
  if (is.null(outname))
    {outname<-"project.mesh.ply"
   }
  if (is.character(tarmesh))
    { target <- tarmesh
    }
  else
    {
      mesh2ply(tarmesh,"target")
    }

  mesh2ply(mesh1,"reference")
  
  if (is.null(mesh1$it))
    {
      cmd <- paste("rayproject reference.ply ",target," -cloud",options," -t ",tol," -o ",outname,sep="")
      system(paste("rayproject reference.ply ",target," -cloud",options," -t ",tol," -o ",outname,sep=""),ignore.stdout=ignore.stdout)
    }
  else
    { cmd <- paste("rayproject reference.ply ",target,options," -t ",tol," -o ",outname,sep="")
      system(paste("rayproject reference.ply ",target,options," -t ",tol," -o ",outname,sep=""),ignore.stdout=ignore.stdout)
    }
  

 # print (cmd)

  outmesh<-ply2mesh(outname,readnormals=TRUE)
  if (clean)
    {unlink(c("reference.ply","target.ply","project.mesh.ply"))
   }
  if (readback)
    {return(outmesh)
   }
}
