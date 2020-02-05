#' visualise differences between two superimposed sets of 3D landmarks
#' 
#' visualise differences between two superimposed sets of 3D landmarks by
#' deforming a cubic grid based on a thin-plate spline interpolation
#' 
#' 
#' @param matrix reference matrix containing 3D landmark coordinates or mesh of class "mesh3d"
#' @param tarmatrix target matrix containing 3D landmark coordinates or mesh of class "mesh3d"
#' @param ngrid number of grid lines to be plotted; ngrid=0 suppresses grid
#' creation.
#' @param align logical: if TRUE, \code{tarmatrix} will be aligned rigidly to \code{matrix}
#' @param lwd width of lines connecting landmarks.
#' @param showaxis integer (vector): which dimensions of the grid to be
#' plotted. Options are combinations of 1,2 and 3.
#' @param show integer (vector): if c(1:2) both configs will be plotted, show = 1 only plots the reference and show = 2 the target
#' @param lines logical: if TRUE, lines between landmarks will be plotted.
#' @param lcol color of lines
#' @param add logical: add to existing rgl window.
#' @param col1 color of "matrix"
#' @param col2 color of "tarmat"
#' @param type "s" renders landmarks as spheres; "p" as points - much faster
#' for very large pointclouds.
#' @param size control size/radius of points/spheres
#' @param pcaxis logical: align grid by shape's principal axes.
#' @param ask logical: if TRUE for > 1000 coordinates the user will be asked to prefer points over spheres.
#' @param margin margin around the bounding box to draw the grid
#' @param createMesh logical: if TRUE, a triangular mesh of spheres and displacement vectors (can take some time depending on number of reference points and grid density).
#' @param slice1 integer or vector of integers: select slice(s) for the dimensions
#' @param slice2 integer or vector of integers: select slice(s) for the dimensions
#' @param slice3 integer or vector of integers: select slice(s) for the dimensions
#' @param gridcol define color of grid
#' @param gridwidth integer: define linewidth of grid
#' @param ... additional parameters passed to \code{\link{rotonto}} in case \code{align=TRUE}
#' @return if \code{createMesh=TRUE}, a mesh containing spheres of reference and target as well as the displacement vectors is returned.
#' @author Stefan Schlager
#' @seealso \code{\link{tps3d}}
#' 
#' @examples
#' \dontrun{
#' data(nose)
#' deformGrid3d(shortnose.lm,longnose.lm,ngrid=10)
#'
#' ## select some slices
#' deformGrid3d(shortnose.lm,longnose.lm,showaxis=1:3,ngrid=10,slice1=2,slice2=5,slice3=7)
#' }
#' @importFrom Rvcg vcgSphere
#' @importFrom rgl translate3d
#' @export
deformGrid3d <- function(matrix,tarmatrix,ngrid=0,align=FALSE,lwd=1,showaxis=c(1, 2), show=c(1,2),lines=TRUE,lcol=1,add=FALSE,col1=2,col2=3,type=c("s","p"), size=NULL, pcaxis=FALSE,ask=TRUE,margin=0.2,createMesh=FALSE,slice1=NULL,slice2=NULL,slice3=NULL,gridcol=1, gridwidth=1,...)
{
    if (inherits(matrix,"mesh3d"))
        matrix <- vert2points(matrix)
    if (inherits(tarmatrix,"mesh3d"))
        tarmatrix <- vert2points(tarmatrix)
    if (align)
        tarmatrix <- rotonto(matrix,tarmatrix,reflection=FALSE,...)$yrot
    type <- type[1]
    if (dim(matrix)[1] > 1000 && type =="s" && (is.null(size) || size > 0) && ask) {
        answer <- readline("You have a lot of landmarks\n Render them as points (faster)? (yes/NO)\n")
        if (! substr(answer,1L,1L) %in% c("n","N"))
            type <- "p"
    }
    out3d <- spheres3d
    if (type == "p")
        out3d <- points3d
    if (!add)
        open3d()
    
    k <- dim(matrix)[1]
    if (is.null(size)) {
        if (type != "p")
            sz <- (cSize(matrix)/sqrt(k)*(1/80))
        else
            sz <- 10
    } else {
        sz <- size
    }
    if (is.numeric(sz) && sz > 0)
        if (1 %in% show)
            out3d(matrix,col=col1,radius=sz, size=sz)
    if(2 %in% show) 
        if (is.numeric(sz) && sz > 0)
            out3d(tarmatrix,col=col2,radius=sz, size=sz)
    if (lines) {
        linemesh <- list()
        linemesh$vb <- t(cbind(rbind(matrix,tarmatrix),1))
        linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
        class(linemesh) <- "mesh3d"
        wire3d(linemesh,lwd=lwd,col=lcol,lit=FALSE)
    }
    x2 <- x1 <- x3 <- c(0:(ngrid-1))/ngrid;
    x0 <- as.matrix(expand.grid(x1,x2,x3))
    
    cent.mat <- scale(matrix, scale=FALSE)
    mean.mat <- colMeans(matrix)
    
    if (ngrid > 1) {
        xrange <- diff(range(matrix[,1]))
        yrange <- diff(range(matrix[,2]))
        zrange <- diff(range(matrix[,3]))
        xrange1 <- diff(range(tarmatrix[,1]))
        yrange1 <- diff(range(tarmatrix[,2]))
        zrange1 <- diff(range(tarmatrix[,3]))
        maxi <- max(c(xrange,yrange,zrange,xrange1,yrange1,zrange1))
        maxi <- (1+margin)*maxi
        x0 <- maxi*x0
        x0 <- scale(x0, scale=FALSE)
        if (pcaxis)
            space <- eigen(crossprod(cent.mat))$vectors
        else
            space <- diag(3)
        x0orig <- t(t(x0%*%space)+mean.mat)
        x0 <- tps3d(x0orig,matrix,tarmatrix,lambda = 1e-8,threads=1)
        
        ## create deformation cube
        outmesh <- list(vb = rbind(t(x0),1))
        class(outmesh) <- "mesh3d"
        
        yinit <- xinit <- zinit <- NULL
        if (2 %in% showaxis) {
            xinit0 <- xinit <- (c(1,2,2+ngrid,1+ngrid))
            for (i in 1:(ngrid-2))
                xinit <- cbind(xinit,(xinit0+i))
            
            xinit0 <- xinit
            for (i in 1:(ngrid-2))
                xinit <- cbind(xinit,xinit0+(i*ngrid))
            
            xinit0 <- xinit
            if (is.null(slice2)) {
                for (i in 1:(ngrid-1))
                    xinit <- cbind(xinit,xinit0+(i*ngrid^2))
            } else {
                xinit <- xinit0+(0:(ngrid-1))[slice2[1]]*ngrid^2
                if (length(slice2) > 1) 
                    for (i in 2:length(slice2))
                        xinit <- cbind(xinit,xinit0+(0:(ngrid-1))[slice2[i]]*ngrid^2)
            }
        }
        if (1 %in% showaxis) {
            yinit0 <- yinit <- c(ngrid,ngrid+ngrid^2, 2*ngrid+ngrid^2, 2*ngrid)
            for( i in 1:(ngrid-2))
                    yinit <- cbind(yinit,yinit0+i*ngrid)
            yinit0 <- yinit
            for (i in 1:(ngrid-2))
                    yinit <- cbind(yinit,yinit0+i*ngrid^2)
            yinit0 <- yinit


            if (is.null(slice1)) {
            for (i in 1:(ngrid-1))
                    yinit <- cbind(yinit,yinit0-i)
            } else {
                yinit <- yinit0-(0:(ngrid-1))[slice1[1]]
                if (length(slice1) > 1) 
                    for (i in 2:length(slice1))
                        yinit <- cbind(yinit,yinit0-(0:(ngrid-1))[slice1[i]])
            }
            
        }
        if (3 %in% showaxis) {
            zinit0 <- zinit <- (c(2,1,1+ngrid^2,2+ngrid^2))
            for( i in 1:(ngrid-2))
                zinit <- cbind(zinit,(zinit0+i))
            
            zinit0 <- zinit
            for (i in 1:(ngrid-2))
                zinit <- cbind(zinit,zinit0+(i*ngrid^2))
            
            zinit0 <- zinit
            if (is.null(slice3)) {
            for (i in 1:(ngrid-1))
                zinit <- cbind(zinit,zinit0+(i*ngrid))
            } else {
                zinit <- zinit0+(0:(ngrid-1))[slice3[1]]*ngrid
                if (length(slice3) > 1) 
                    for (i in 2:length(slice3))
                        zinit <- cbind(zinit,zinit0+(0:(ngrid-1))[slice3[i]]*ngrid)
            }
        }
        outmesh$ib <- cbind(xinit,yinit,zinit)
        wire3d(outmesh,lit=F,col=gridcol,lwd=gridwidth)
    }
    ## create a mesh displaying the stuff
    if (createMesh) {
        mysphere <- vcgSphere(subdivision=1)
        mysphere <- scalemesh(mysphere,sz,"none")
        col1mesh <- rgb(t(col2rgb(col1)), maxColorValue = 255)
        matmesh <- lapply(1:nrow(matrix),function(x) x <- mysphere)
        matmesh <- lapply(1:nrow(matrix),function(x) x <- translate3d(matmesh[[x]],x=matrix[x,1],y=matrix[x,2],z=matrix[x,3]))
        matmesh <- mergeMeshes(matmesh)
        matmesh$material$color <- rep(col1mesh,ncol(matmesh$vb))

        col2mesh <- rgb(t(col2rgb(col2)), maxColorValue = 255)
        tarmatmesh <- lapply(1:nrow(tarmatrix),function(x) x <- mysphere)
        tarmatmesh <- lapply(1:nrow(tarmatrix),function(x) x <- translate3d(tarmatmesh[[x]],x=tarmatrix[x,1],y=tarmatrix[x,2],z=tarmatrix[x,3]))
        tarmatmesh <- mergeMeshes(tarmatmesh)
        tarmatmesh$material$color <- rep(col2mesh,ncol(matmesh$vb))

        allMerge <- mergeMeshes(matmesh,tarmatmesh)
        if (lines) {
            lcolmesh <- rgb(t(col2rgb(lcol)), maxColorValue = 255)
            diffs <- tarmatrix-matrix
            dists <- sqrt(rowSums(diffs^2))
            mylinemesh <- mergeMeshes(lapply(1:nrow(matrix),function(x) cylinder(matrix[x,],diffs[x,],dists[x],fine=10,radius=lwd*(sz/10))))
            mylinemesh$material$color <- rep(lcolmesh,ncol(mylinemesh$vb))
            allMerge <- mergeMeshes(allMerge,mylinemesh)
        }
        if (ngrid > 0) {
            cageverts <- vert2points(outmesh)
            refinds <- as.vector(t(outmesh$ib))
            tarinds <- as.vector(t(outmesh$ib[c(2:4,1),]))
            cagetarmatrix <- cageverts[tarinds,]
            cagematrix <- cageverts[refinds,]
            cagediffs <- cageverts[tarinds,]-cageverts[refinds,]
            cagedists <- sqrt(rowSums(cagediffs^2))
            mycagemesh <- mergeMeshes(lapply(1:nrow(cagematrix),function(x) cylinder(cagematrix[x,],cagediffs[x,],cagedists[x],fine=4,radius=lwd*(sz/10))))
            mycagemesh$material$color <- rep("#000000",ncol(mycagemesh$vb))
            #mycagemesh <- tps3d(mycagemesh,matrix,tarmatrix,lambda = 1e-8,threads=1)
            allMerge <- mergeMeshes(allMerge,mycagemesh)
        }
        invisible(allMerge)
    }
    
}


cylinder <- function(x,dirs,length,radius=1,fine=20,adNormals=FALSE)
    {
### create a 3D mesh representing a cylinder and place it in a requested positon
### create initial circle ###
        seqby <- 2/fine
        circ <- seq(from=0,to=(2-seqby), by=seqby)*pi
        data <- matrix(0,fine,3)
        data[,1] <- sin(circ)
        data[,2] <- cos(circ)
        data <- data*radius
        lc <- dim(data)[1]
        
###  begin create faces of a cylinder ###
        it <-NULL
        for (i in 1:(lc-1)) {
            face0 <- c(i,i+1,i+lc)
            face1 <- c(i+1,i+lc+1,i+lc)
            it <- rbind(it,face0,face1)
        }
        it <- rbind(it,c(lc,1,lc+1))
        it <- rbind(it,c(2*lc,lc,lc+1))
### close lids ###
        for (i in 1:(lc-1))
            it <- rbind(it,c(lc*2+1,i+1,i))
        it <- rbind(it,c(1,lc,lc*2+1))
        
        for (i in (lc+1):(2*lc-1))
            it <- rbind(it,c(lc*2+2,i,i+1))
        it <- rbind(it,c(lc*2+2,lc*2,lc+1))
### end faces creation ###
        
### rotate initial circle and create opposing circle ###    
        dirs <- dirs/sqrt(sum(dirs^2))
        yz <- tangentPlane(dirs)
        rmat <- cbind(yz$z,yz$y,dirs)
        data <- t(rmat%*%t(data))
        data2 <- t(apply(data,1,function(x){x <- x+length*dirs;return(x)}))

### create cylinder mesh ###
        datavb <- rbind(data,data2)
        x0 <- c(0,0,0)
        x1 <- x0+length*dirs
        datavb <- rbind(datavb,x0,x1)
        cylvb <-  t(cbind(datavb,1))
        colnames(cylvb) <- NULL
        cyl <- list()
        class(cyl) <- "mesh3d"
        cyl$vb <- cylvb
        cyl$it <- t(it)
        cyl <- invertFaces(cyl)
        cyl$normals <- NULL
        cyl <- translate3d(cyl,x[1],x[2],x[3])
        if (adNormals)
            cyl <- vcgUpdateNormals(cyl)
            
        return(cyl)
    }
