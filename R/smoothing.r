smoothMesh3d <- function(mesh,method=c("taubin","laplace","HClaplace"),iteration=10)
{
  option <- NULL
  method <- substring(method[1],1L,1L)
  mesh2ply(mesh,"smoothdump")
  if (method == "l")
    {
      option <- "--laplace "
    }
  else if (method %in% c("H","h"))
   {
      option <- "--lapHC "
    }
    command <- paste("trismooth smoothdump.ply -it ",iteration," ",option,"smoothdump.ply",sep="")
  system(command)
  mesh <- ply2mesh("smoothdump.ply")
  unlink("smoothdump.ply")
  invisible(mesh)
}
