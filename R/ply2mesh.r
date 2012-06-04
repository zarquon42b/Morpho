ply2mesh<-function (filename, adnormals = TRUE,readnormals=FALSE,readcol=FALSE)
{

  x <- filename
  A <- readLines(x, n = 100)
  end <- which(A == "end_header")
  infos <- A[1:end]
  vertinfo <- strsplit(A[grep("element vertex", infos)], " ")
  faceinfo <- strsplit(A[grep("element face", infos)], " ")
  texinfo<-NULL
  colmat <- NULL
  material <- NULL
                                        #if (length(grep("property list uchar float texcoord",A))==1) 
  qualinfo<-grep("property float quality", infos)
  vertbegin <- grep("element vertex",infos)
  facebegin <- grep("element face",infos)
  if (length(qualinfo)==1)
    {qualine<-qualinfo-vertbegin
     qual<-TRUE
   }
  else 
    {qual <- FALSE}
  fn <- as.numeric(faceinfo[[1]][3])
  vn <- as.numeric(vertinfo[[1]][3])
  vert.all <- read.table(x, skip = end, sep = " ", nrows = vn)
  vert <- apply(vert.all[, 1:3], 2, as.numeric)
  vert.n <- NULL
  quality<-NULL
  if (length(grep("property float nx", infos)) == 1)
    {
      normstart <- grep("property float nx", infos)-vertbegin
      vert.n <- t(vert.all[, normstart:(normstart+2)])
      if (qual)
        {quality<-as.vector(vert.all[,qualine])}
    }

  if (readcol)##check for colored vertices
    {
      color <- grep("property uchar red",infos)
      if (length(color > 0))
          { if (color[1] < facebegin)
              {
                colbegin <- color[1]-vertbegin
                colmat <- vert.all[,colbegin:(colbegin+2)]
                colmat <- rgb(colmat[,1],colmat[,2],colmat[,3],maxColorValue=255)
              }
          }
    }

  if (fn !=0)
    {
      face.all <- read.table(x, skip = end + vn, nrows = fn)
      face <- t(face.all[, 2:4]+1)
      
      if (!is.null(colmat))
        {
          dimface <- dim(face)
          colfun <- function(x)
            {
              x <- colmat[x]
              return(x)
            }
          material$color <- colfun(face)
        }
      mesh <- list(vb = rbind(t(vert), 1), it = face, primitivetype = "triangle", material = material,normals = vert.n)      
      class(mesh) <- c("mesh3d", "shape3d")
    }
  else
    {if (is.null(vert.n) || readnormals==FALSE)
       {cat(paste("mesh contains no faces. Vertices will be stored in a",vn,"x 3 matrix\n"))
        mesh<-vert}	
    else if (readnormals)
      {cat(paste("mesh contains no faces. vertices and vertex normals are stored in a list\n"))
       mesh<-list(vb=t(vert),normals=vert.n)
     }	
     
   }
### generate object of class mesh3d ###	
  
### add TexCoords ###
  if (length(grep("property list uchar float texcoord",A))==1 && length(grep("comment TextureFile",A))==1)
    {texn<-face.all[1,5]
     tex<-face.all[,c(6:(6+(texn-1)))]
     mesh$tex<-t(tex)
     mesh$TextureFile<-strsplit(A[grep("comment TextureFile", infos)], " ")[[1]][3]
   } 
    
### check for normals and update if required ###
  if (fn !=0)	
    {if (adnormals && is.null(mesh$normals) ) {
      cat("calculating normals...\n")
      mesh <- adnormals(mesh)	
    }
   }
  if (qual && readnormals)
    {mesh$quality<-quality
   }
  return(mesh)
}

