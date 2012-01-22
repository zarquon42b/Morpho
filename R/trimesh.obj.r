trimesh.obj<-function(obj,wire="s",col=1,size=size,plot=T,shapeout=T,addNormals=TRUE)     #options: w=wireframe, s=shade, p=points
{     	
	
	vert<-obj[which(obj[,1]=="v"),1:4]
        face<-obj[which(obj[,1]=="f"),1:4]
      	#vn<-obj[which(obj[,1]=="vn"),1:4]
	face.mat<-as.matrix(face[,2:4])
	vert.mat<-as.matrix(vert[,2:4])
	#vert.ind<-which(obj[,1]=="v")
	#vn.ind<-which(obj[,1]=="vn")	
	
	mesh<-tmesh3d(t(vert.mat),t(face.mat),homogeneous=FALSE)
		if (addNormals && is.null(mesh$normal))
			{
			mesh<-addNormals(mesh)
			}

        if (wire=="w"&& plot==TRUE)
            {wire3d(mesh,color=col)}

        if (wire=="s" && plot==TRUE)
            {shade3d(mesh,color=col)}
        if  (wire=="p" && plot==TRUE)
            {dot3d(mesh,color=col)}


        if (shapeout==T){return(mesh=mesh)}

}
