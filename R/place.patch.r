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
plotAtlas <- function(atlas, pt.size=NULL, alpha=1, render=c("w","s"), point=c("s", "p"), meshcol="white", add=TRUE, legend=TRUE)
    {
        outid <- NULL
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        point <- point[1]
        ## set point/sphere sizes
        radius <- pt.size
        if (is.null(radius)) {
            if (point == "s")
                radius <- (cSize(atlas$landmarks)/sqrt(nrow(atlas$landmarks)))*(1/30)
            else
                radius <- 10
        }
        size <- radius
        render <- render[1]
        if (point == "s") {
            rendpoint <- spheres3d
        } else if (point == "p") {
            rendpoint <- points3d
        } else {
            stop("argument \"point\" must be \"s\" for spheres or \"p\" for points")
        }
        if (render=="w") {
            rend <- wire3d
        } else {
            rend <- shade3d
        }
        if (add)
            open3d()
        if (!is.null(atlas$mesh))
            outid <- rend(atlas$mesh, col=meshcol, alpha=alpha)
        outid <- c(outid, rendpoint(atlas$landmarks,col=2, radius=radius, size=size))
        outid <- c(outid,rendpoint(atlas$patch,col=3,radius=radius/2, size=size/2))
        if (!is.null(atlas$corrCurves))
            outid <- c(outid, rendpoint(atlas$landmarks[unlist(atlas$corrCurves),],col=4,radius=radius+0.001, size=size+1))
        if (!is.null(atlas$patchCurves))
            outid <- c(outid,rendpoint(atlas$patch[unlist(atlas$patchOutlines),],col=5,radius=radius/2+0.001,size=(size/2)+1))
        if (legend) {
            plot(0,0, xlab="", ylab="", axes =F, cex=0,xlim=c(-1,1), ylim=c(-1,1))
            legend(-1,1, pch=20, cex=2, col=2:5, legend=c("landmarks", "patch", "curves on all specimen", "curves only on atlas"))
        }
        invisible(outid)
    }

placePatch <- function(atlas, dat.array, path, prefix=NULL, fileext=".ply", ray=TRUE, inflate=NULL,tol=inflate, relax.patch=TRUE, keep.fix=NULL, rhotol=NULL, silent=FALSE)
    {
        if (!inherits(atlas, "atlas"))
            stop("please provide object of class atlas")
        if (is.null(keep.fix))
            keep.fix <- 1:dim(atlas$landmarks)[1]
        if (is.null(tol) && !is.null(inflate))
            tol <- inflate
        
        patched <- place.patch(dat.array, path, atlas.mesh =atlas$mesh, atlas.lm = atlas$landmarks, patch =atlas$patch, curves=atlas$patchCurves, prefix=prefix, tol=tol, ray=ray, outlines=atlas$corrCurves, inflate=inflate, relax.patch=relax.patch, rhotol=rhotol, fileext=fileext, SMvector = keep.fix, silent=silent)
        return(patched)
    }

place.patch <- function(dat.array,path,atlas.mesh,atlas.lm,patch,curves=NULL,prefix=NULL,tol=5,ray=T,outlines=NULL,SMvector=NULL,inflate=NULL,relax.patch=TRUE,rhotol=NULL,fileext=".ply", silent=FALSE)
    {
        k <- dim(dat.array)[1]
        deselect=TRUE
        fix <- which(c(1:k) %in% SMvector)
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
                tmp.data <- projRead(dat.array[,,i],tmp.name,readnormals=TRUE, ignore.stdout=silent)
            else
                tmp.data <- projRead(dat.array,tmp.name,readnormals=TRUE, ignore.stdout=silent)
            
### relax existing curves against atlas ###
            if (!is.null(outlines)) {
                sm <- SMvector
                U <- .calcTang_U_s(t(tmp.data$vb[1:3,]),t(tmp.data$normals[1:3,]),SMvector=SMvector,outlines=outlines,surface=NULL,deselect=deselect)
                slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
                tmp.data <- projRead(slide,tmp.name,readnormals=TRUE, ignore.stdout=silent)
                tps.lm <- tps3d(patch,atlas.lm,slide)
            } else if (!is.null(SMvector) && is.null(outlines)) {
                sm <- SMvector
                tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            } else {
                sm <- 1:k
                tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            }

            slide <- t(tmp.data$vb[1:3,])
            slidenormals <- t(tmp.data$normals[1:3,])
            slide[fix,] <- dat.array[fix,,i] #replace projected points with original for fix landmarks
### use for mullitlayer meshes to avoid projection inside
            if (!is.null(inflate)) {
                atlas.warp <- warp.mesh(atlas.mesh,atlas.lm,slide, silent=silent)
                tps.lm <- projRead(tps.lm,atlas.warp,readnormals=TRUE,smooth=TRUE, ignore.stdout = silent)
                warp.norm <- tps.lm$normals[1:3,]### keep projected normals
                
                tps.lm$vb[1:3,] <- tps.lm$vb[1:3,]+inflate*tps.lm$normals[1:3,] ###inflate outward along normals
                tps.lm <- ray2mesh(tps.lm,tmp.name,inbound=TRUE,tol=tol,angmax=rhotol, ignore.stdout=silent) ### deflate in opposite direction
            } else {## just project warped patch on surface (suitable for singlelayer meshes)
                tps.lm <- projRead(tps.lm,tmp.name,readnormals=TRUE, ignore.stdout=silent )
            }
            
            relax <- rbind(slide,t(tps.lm$vb[1:3,]))
            normals <- rbind(slidenormals,t(tps.lm$normals[1:3,]))
            
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
                
                U1 <- .calcTang_U_s(relax, normals,SMvector=sm,outlines=outltmp,surface=surface,free=free,deselect=deselect)
                tps.lm <- calcGamma(U1$Gamma0,L1$Lsubk3,U1$U,dims=3)$Gamatrix[c((k+1):(patch.dim+k)),]
                tps.lm <- projRead(tps.lm,tmp.name,readnormals=FALSE, ignore.stdout = silent)
            } else {# end relaxation ########################
                tps.lm <- t(tps.lm$vb[1:3,])
            }
            
            if (!usematrix)
                out[,,i] <- rbind(dat.array[,,i],tps.lm)
            else
                out <- rbind(dat.array,tps.lm)
        }
        return(out)
    }

