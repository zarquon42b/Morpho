rotmesh.onto <- function(mesh, refmat, tarmat, adnormals=FALSE, scale=FALSE, reflection=FALSE)
{
  rot <- rotonto(tarmat,refmat,scale=scale,reflection=reflection)
  mesh$vb[1:3,] <- mesh$vb[1:3,]-rot$transy
  mesh$vb[1:3,] <- t(t(mesh$vb[1:3,])%*%rot$gamm)
  if (scale)
      mesh$vb[1:3,] <- mesh$vb[1:3,]*rot$bet
    
  mesh$vb[1:3,] <- mesh$vb[1:3,]+rot$trans
  if (sign(det(rot$gamm) < 0 && reflection))
      mesh <- conv2backf(mesh)
  if (adnormals) 
      mesh <- adnormals(mesh)
  if (!is.null(mesh$normals) && !adnormals)
      mesh$normals[1:3,] <- t(t(mesh$normals[1:3,]) %*% rot$gamm)
   
  return(list(mesh=mesh,yrot=rot$yrot))
}
