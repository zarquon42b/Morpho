mesh2ply<-function(x,filename="default",col=NULL,writeNormals=FALSE)
{	

  if (is.matrix(x))
    {x<-list(vb=x)}
  filename<-paste(filename,".ply",sep="")
  vert<-x$vb[1:3,]
  vert<-round(vert,digits=6)
  if (!is.null(x$it))
    {face<-x$it-1
     fd<-3
     fn<-dim(face)[2]
   }
  else 
    {fn<-0
   }	
  vert.all<-vert
  vn<-dim(vert)[2]
  vn.all<-3
  texfile<-x$TextureFile

  if (is.null(col) && !is.null(x$material$color))
    {
      col=rep("#FFFFFF",vn)
      tmp1 <- data.frame(it=as.vector(x$it))
      tmp1$rgb <- as.vector(x$material$color)
      tmp1 <- unique(tmp1)
      col[tmp1$it] <- tmp1$rgb
    }

  if (!writeNormals)
    {
      x$normals <- NULL
    }
  
### start writing to file ###
  
### write header ###
  cat("ply\nformat ascii 1.0\ncomment MORPHO generated\n",file=filename)
  
### check for Texture information and write to header ###	
  
  if (is.character(texfile))
    {cat(paste("comment TextureFile ",texfile,"\n",sep=""),file=filename,append=TRUE) 
   }
  
### write vertex infos to header ###
  
	v.info<-paste("element vertex ",vn,"\n","property float x\nproperty float y\nproperty float z\n",sep="")
 
	cat(v.info,file=filename,append=TRUE)	
		if (!is.null(x$normals))
			{cat("property float nx\nproperty float ny\nproperty float nz\n",file=filename,append=TRUE)
			norma<-round(x$normals[1:3,],digits=6)
			vert.all<-rbind(vert,norma)	
			vn.all<-6		
			}
   if (!is.null(col))    
     {
      
       v.info <- cat("property uchar red\nproperty uchar green\nproperty uchar blue\n",file=filename,append=T)
     }
	### write face infos and texture infos to header and finish ###	
	
	cat(paste("element face ",fn,"\n",sep=""),file=filename,append=TRUE)	
	if(!is.null(x$tex) && is.character(texfile))
		{cat("property list uchar int vertex_indices\nproperty list uchar float texcoord\nend_header\n",file=filename,append=TRUE)	
		}
	else
		{cat("property list uchar int vertex_indices\nend_header\n",file=filename,append=TRUE)	
		}
	
	### write vertices and vertex normals ###
  vert.all <- data.frame(t(vert.all))
  if (!is.null(col))
    {
      if (is.vector(col))
        {
          colout <- t(col2rgb(col))
        }
      else
        {
          colout <- matrix(col2rgb(col),vn,3,byrow=T)
        }
      vert.all <- cbind(vert.all,colout)
    }
  
      write.table(format(vert.all,scientific=F,trim=T),file=filename,sep=" ",append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")	

	
	### write face and Texture information ###	
	if (!is.null(x$it)){
	if(!is.null(x$tex) && !is.null(x$TextureFile))
		{tex<-t(x$tex)
		texn<-dim(tex)[2]
		faceraw<-rbind(fd,face)
		facetex<-t(cbind(texn,tex))
		write.table(format(t(rbind(faceraw,facetex)),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")			
		}
	else
		{write.table(format(t(rbind(fd,face)),scientific=F,trim=T),file=filename,sep=" ",append=TRUE,quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")	
		}
	}
}
