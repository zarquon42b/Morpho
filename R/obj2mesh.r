obj2mesh <- function(filename,adnormals=TRUE)
{	
	obj <- read.obj(filename)
	vert <- obj[which(obj[,1]=="v"),1:4]
	        
	face <- obj[which(obj[,1]=="f"),1:4]
      	vn <- obj[which(obj[,1]=="vn"),1:4]
	face.mat <- as.matrix(face[,2:4])
	vert.mat <- apply(vert[,2:4],2,as.numeric)
	#vert.ind_old <- which(obj[,1]=="v") 	# original index of vertices
	#vn.ind <- which(obj[,1]=="vn")  		# original index of vertex normals
	
	#vn.mat <- apply(vn[,2:4],2,as.numeric)
	#vn.mat <- as.matrix(vn.mat)		
	#vn.mat_new <- vert.mat	
	
	if (length(grep("//",face.mat[1,1]))!=0)
		{	
		
       		write.table(face.mat,file="facedump", quote = F, row.names = FALSE, col.names = FALSE, na = "",sep="//")
		face.mat <- read.table("facedump",sep="/")[,c(1,5,9)]	
		unlink("facedump")	
		}
	
	
		
	#if (dim(vn)[1]!=0 && sum(abs(vn.mat)) > 0 ) ### check for valid vertex normals
	#	{for (i in 1:dim(vert.mat)[1])
	#		{ptr <- vn.ind[max(which(vn.ind < vert.ind_old[i]))]
	#		vn.mat_new[i,] <- vn.mat[ptr,]
	#		}
	#	vn.mat_new <- t(vn.mat_new)
	#	}
	#else 
	#	{vn.mat_new <- NULL}
        mesh <- list()
        class(mesh) <- "mesh3d"
	mesh$vb <- rbind(t(vert.mat),1)
        mesh$it <- t(face.mat)
	#mesh$normals <- vn.mat_new
	
	if (adnormals && is.null(mesh$normals))
			{
			mesh <- adnormals(mesh)
			}
	return(mesh)
}
