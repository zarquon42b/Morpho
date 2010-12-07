bone2skin<-function(data,bonemesh,skinmesh,meshtype=".ply")
{		out<-list()
		out$vb<-t(data)
		proj.back(data,bonemesh)
		system(paste("trinorm_project out_cloud.ply ", skinmesh,sep=""))
		out<-ply2mesh("out_norm.ply")
		return(out)
}
		
		
