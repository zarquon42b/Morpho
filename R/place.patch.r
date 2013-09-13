createAtlas <- function(mesh, landmarks, patch, corrCuves=NULL, patchCurves=NULL)
    {
        atlas <- list()
        class(atlas) <- "atlas"
        atlas$mesh <- mesh
        atlas$landmarks <- landmarks
        atlas$patch <- patch
        atlas$corrCurves <- corrCuves
        atlas$patchCurves<- patchCurves
        return(atlas)
    }
plotAtlas <- function(atlas, radius=1, meshcol="white", add=TRUE)
    {
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        if (add)
            open3d()
        wire3d(atlas$mesh, col=meshcol)
        spheres3d(atlas$landmarks,col=2)
        spheres3d(atlas$patch,col=3,radius=radius/2)
        if (!is.null(atlas$corrCurves))
            spheres3d(atlas$landmarks[unlist(atlas$corrCurves),],col=4,radius=radius+0.001)
        if (!is.null(atlas$patchCurves))
            spheres3d(atlas$patch[unlist(atlas$patchOutlines),],col=5,radius=radius/2+0.001)
        plot(0,0, xlab="", ylab="", axes =F, cex=0,xlim=c(-1,1), ylim=c(-1,1))
        legend(-1,1, pch=20, cex=2, col=2:5, legend=c("landmarks", "patch", "curves on all specimen", "curves only on atlas"))
    }

placePatch <- function(atlas, dat.array, path, prefix=NULL, fileext=".ply", tol=5, ray=TRUE, inflate=NULL, relax.patch=TRUE, keep.fix=NULL, rhotol=NULL)
    {
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        if (is.null(keep.fix))
            keep.fix <- 1:dim(atlas$landmarks)[1]
            
        patched <- place.patch(dat.array, path, atlas.mesh =atlas$mesh, atlas.lm = atlas$landmarks, patch =atlas$patch, curves=atlas$patchCurves, prefix=prefix, tol=tol, ray=ray, outlines=atlas$corrCurves, inflate=inflate, relax.patch=relax.patch, rhotol=rhotol, fileext=fileext,SMvector = keep.fix)
        return(patched)
    }

place.patch <- function(dat.array,path,atlas.mesh,atlas.lm,patch,curves=NULL,prefix=NULL,tol=5,ray=T,outlines=NULL,SMvector=NULL,deselect=TRUE,inflate=NULL,relax.patch=TRUE,rhotol=NULL,fileext=".ply")
    {
        k <- dim(dat.array)[1]
        patch.dim <- dim(patch)[1]
        usematrix <- FALSE
        if (! is.matrix(dat.array)) {
            n <- dim(dat.array)[3]
            out <- array(NA,dim=c((patch.dim+k),3,n))
            dimnames(out)[[3]] <- dimnames(dat.array)[[3]]
            name <- dimnames(dat.array)[[3]]
        } else {
            usematrix <- TRUE
            n <- 1
            out <- matrix(NA,(patch.dim+k),3)
            name <- NULL
        }
        
        L <- CreateL(atlas.lm)
        L1 <- CreateL(rbind(atlas.lm,patch))
        
        for(i in 1:n) {
            tmp.name <- paste(path,prefix,name[i],fileext,sep="")
            if (!usematrix)
                tmp.data <- projRead(dat.array[,,i],tmp.name,readnormals=TRUE)
            else
                tmp.data <- projRead(dat.array,tmp.name,readnormals=TRUE)
            
### relax existing curves against atlas ###
            if (!is.null(outlines)) {
                sm <- SMvector
                U <- .calcTang_U_s(t(tmp.data$vb[1:3,]),t(tmp.data$normals[1:3,]),SMvector=SMvector,outlines=outlines,surface=NULL,deselect=deselect)
                slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
                slide <- projRead(slide,tmp.name,readnormals=FALSE)
                tps.lm <- tps3d(patch,atlas.lm,slide)
            } else if (!is.null(SMvector) && is.null(outlines)) {
                sm <- SMvector
                slide <- t(tmp.data$vb[1:3,])           
                tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            } else {
                sm <- 1:k
                slide <- t(tmp.data$vb[1:3,])           
                tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            }
### use for mullitlayer meshes to avoid projection inside
            if (!is.null(inflate)) {
                atlas.warp <- warp.mesh(atlas.mesh,atlas.lm,slide)
                tps.lm <- projRead(tps.lm,atlas.warp,readnormals=TRUE,smooth=TRUE)
                warp.norm <- tps.lm$normals[1:3,]### keep projected normals
                
                tps.lm$vb[1:3,] <- tps.lm$vb[1:3,]+inflate*tps.lm$normals[1:3,] ###inflate outward along normals
                tps.lm <- ray2mesh(tps.lm,tmp.name,inbound=TRUE,tol=tol,angmax=rhotol) ### deflate in opposite direction
                relax <- rbind(slide,t(tps.lm$vb[1:3,]))
                normals <- rbind(slide,t(tps.lm$normals[1:3,]))
                
                surface <- c((k+1):(patch.dim+k))  ## define surface as appended to preset landmarks
                free <- NULL
### compare normals of projection and original points
                if (!is.null(rhotol)) {
                    rho <- NULL
                    for (j in 1:patch.dim)
                        rho[j] <- angle.calc(tps.lm$normals[1:3,j],warp.norm[1:3,j])
                    
                    rhoex <- which(rho > rhotol) 
                    if (length(rhoex) > 0) {
                        free <- surface[rhoex]
                        surface <- surface[-rhoex]
                    }
                }
                gc()
### end compare normals #### 
                
### relax patch against reference ###
                if (relax.patch){ ### relax against reference
                    outltmp <- append(outlines,curves) ## add curves from patch to predefined curves
                    remout <- which(surface %in% curves)
                    
                    if (length(remout) > 0)
                        surface <- surface[-remout] ### remove patch curves from surface 
                    if (length(surface)==0)
                        surface <- NULL
                    
                    U1 <- .calcTang_U_s(relax,normals,SMvector=sm,outlines=outltmp,surface=surface,free=free,deselect=deselect)
                    tps.lm <- calcGamma(U1$Gamma0,L1$Lsubk3,U1$U,dims=3)$Gamatrix[c((k+1):(patch.dim+k)),]
                    tps.lm <- projRead(tps.lm,tmp.name,readnormals=FALSE)
                } else {# end relaxation ########################
                    tps.lm <- t(tps.lm$vb[1:3,])
                }
            } else {## just project warped patch on surface (suitable for singlelayer meshes)
                tps.lm <- projRead(tps.lm,tmp.name,readnormals=FALSE)
            }
            if (!usematrix)
                out[,,i] <- rbind(dat.array[,,i],tps.lm)
            else
                out <- rbind(dat.array,tps.lm)
        }
        return(out)
    }

