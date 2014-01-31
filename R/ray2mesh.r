#' @rdname mesh2mesh
#' @param targetdump specify name of dumped target mesh
#' @param refdump specify name of dumped reference mesh
#' @export
ray2mesh <- function(mesh1,tarmesh,tol=1,angmax=NULL,clean=TRUE,outname=NULL,readback=TRUE,inbound=FALSE,strict=FALSE,ignore.stdout=FALSE,mindist=FALSE,targetdump="target",refdump="reference")
{ 

  options <- NULL
  opt <- FALSE
  target <- paste0(targetdump, ".ply")
  reference <- paste0(refdump,".ply")
  if (inbound == TRUE) {
      opt <- TRUE
      options <- "--inbound" # set option to search along negative normals first
  }
  if (strict == TRUE) {
      opt <- TRUE
      options <- paste(options,"--strict") #mark vertices that are not hit along rays
  }
  if (mindist == TRUE) {
      opt <- TRUE
      options <- paste(options,"--minray") #use closest point in and outward
  }
  if (!is.null(angmax)) {
      opt <- TRUE
      options <- paste(options,"--angmax",angmax) #check for normals
  }
  if (opt) 
      options <- paste(" ",options,sep="")
  if (is.null(outname)) 
      outname <- "project.mesh.ply"
  
  if (is.character(tarmesh)) {
      target <- tarmesh
  } else {
      mesh2ply(tarmesh,targetdump)
  }
  if (is.character(mesh1))
      reference <- mesh1
  else
      mesh2ply(mesh1,refdump)
  
  mesh2ply(mesh1,refdump)
  
  if (is.null(mesh1$it)) {
      #cmd <- paste("rayproject ",reference," ",target," -cloud",options," -t ",tol," -o ",outname,sep="")
      system(paste("rayproject ",reference," ",target," -cloud",options," -t ",tol," -o ",outname,sep=""),ignore.stdout=ignore.stdout)
  } else {
      #cmd <- paste("rayproject reference.ply ",target,options," -t ",tol," -o ",outname,sep="")
      system(paste("rayproject ", reference," ",target," ",options," -t ",tol," -o ",outname,sep=""),ignore.stdout=ignore.stdout)
  }
  
  outmesh <- ply2mesh(outname,readnormals=TRUE, silent=ignore.stdout)
  if (clean) {
      unlink(c(reference,outname))
      if (inherits(tarmesh,"mesh3d"))
          unlink(target)
  }
  if (readback)
      return(outmesh)
}
