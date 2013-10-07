projBack <- function(data,surface,dataname=NULL,outname=NULL,smooth=TRUE,ignore.stdout=FALSE,sign=FALSE)
{
  checkCLI <- try(system("trimesh_project",intern=TRUE),silent=TRUE)
   if(class(checkCLI) == "try-error")
      stop("please install trimesh-tools")
  options <- NULL
  
  if (!smooth)
      options <- paste(options,"--nosmooth")
  if (sign)
      options <- paste(options,"--sign")
  if (is.null(dataname))
      dataname <- "out"
  write.obj(cbind("v",data),filename=dataname)
  if (is.null(outname))
      command <- paste("trimesh_project"," ",dataname,".obj"," ",surface,options,sep="")  
  else
      command <- paste("trimesh_project"," ",dataname,".obj"," ",surface," -o ",outname,options,sep="")
  
  system(command,ignore.stdout=ignore.stdout)
  unlink(paste(dataname,".obj",sep="")) #clean up
}

projRead <- function(lm,mesh,readnormals=TRUE,clean=TRUE,smooth=TRUE,ignore.stdout=FALSE,sign=FALSE,lmdump=NULL,prodump=NULL)
{

  checkCLI <- try(system("trimesh_project",intern=TRUE),silent=TRUE)
   if(class(checkCLI) == "try-error")
      stop("please install trimesh-tools")
  if (is.null(prodump))
    prodump <- "out_cloud.ply"
  if (is.null(lmdump))
    lmdump <- "out"
  if (is.character(mesh))
      projBack(lm,mesh,ignore.stdout=ignore.stdout,sign=sign,smooth=smooth,outname=prodump,dataname=lmdump)
  else {
      dumpname <- paste(prodump,"dump0",sep="")
      dumpfile <- paste(dumpname,".ply",sep="")
      mesh2ply(mesh,dumpname)
      projBack(lm,dumpfile,smooth=smooth,ignore.stdout=ignore.stdout,sign=sign,outname=prodump,dataname=lmdump)
      unlink(dumpfile)
  }
  
  data <- ply2mesh(prodump,readnormals=readnormals)
  if (clean)
      unlink(prodump)
  
  return(data)
}
