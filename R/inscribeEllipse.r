check_inner <- function(poly_x,poly_y,a,b,x0,y0 ){
    minrr = 5
    for (i in (1:length(poly_x))) {
        x=poly_x[i]
        y=poly_y[i]        
        rr  <- ((x-x0)/a)^2+((y-y0)/b)^2
        if (rr < minrr)
            minrr = rr
                                        #if rr < 1 : 
                                        #    badpoint.append([x,y])
    }
    return (minrr)
    
}

get_jumperpoint <- function(poly_x,poly_y,sizeup,a,b,x0,y0   ) {
    badx = NULL
    bady = NULL
    for (i in 1:length(poly_x)) {
        x=poly_x[i]
        y=poly_y[i]        
        rr = ((x-x0)/(a+sizeup))^2+((y-y0)/(b+sizeup))^2
        if (rr < 1) {
            badx <- c(badx,x)
            bady <- c(bady,y)
        }
    }
    xmean = mean(badx)
    ymean = mean(bady)   
    xmove = x0-xmean
    ymove = y0-ymean
    return(c(xmove,ymove))
}



inscribeEllipseOld <- function(poly,step=0.3,iters=90) {
  
    px_old = poly[,1]
    py_old = poly[,2]
    init_point = colMeans(poly)
    init_radius = step


    px = NULL
    py = NULL
    for (i in (1:(length(px_old)-1))) {
        dx = px_old[i+1]-px_old[i]
        dy = py_old[i+1]-py_old[i]
        len1 = (dx*dx+dy*dy)^0.5
        px <- c(px,px_old[i])
        py <- c(py,py_old[i])
        if (len1 >= step){
            count = as.integer(len1/step)
            for (ii in (1:count)) {
                ##print px_old[i]+ 1.0*ii/count*dx
                ##print py_old[i]+ 1.0*ii/count*dy
                px <- c(px,(px_old[i]+ 1.0*ii/count*dx))
                py <- c(py, (py_old[i]+ 1.0*ii/count*dy))
            }
        }
    }
    
    px <- c(px,(px_old[-1]))
    py <- c(py,(py_old[-1]))
    minrr = 6

    rx = init_radius
    ry = init_radius
    xc = init_point[1]
    yc = init_point[2]

    maxarea = step*step
for (iterat in (1:iters)) {
 ##   print ("iterate start %s (%s,%s)->(%s,%s)" % (iterat,xc,yc,rx,ry))
    s1 = step
    s2 = step
    s3 = step
    #draw_ellips(ax,xc,yc,rx,ry,'blue',0.5)

    rxprev = rx 
    ryprev = ry
    xcprev = xc
    ycprev = yc
    if (check_inner(px,py,rx+s1,ry+s2,xc,yc) >= 1) { 
        rx = rx+s1
        ry = ry+s2
    }
    # enlarge
    if (check_inner(px,py,rx+s1,ry,xc,yc) >= 1) 
        rx = rx+s1

    if (check_inner(px,py,rx,ry+s1,xc,yc) >= 1)
        ry = ry+s1

    #distortion
    if (check_inner(px,py,rx+s1,ry-s2,xc,yc) >= 1 && (rx+s1)*(ry-s2) > rx*ry) {
        rx = rx+s1
        ry = ry-s2      
    }
    if (check_inner(px,py,rx-s1,ry+s2,xc,yc) >= 1 && (rx-s1)*(ry+s2) > rx*ry) {
        rx = rx-s1
        ry = ry+s2
    }
    #shift @ enlarge (right dir )
    if (check_inner(px,py,rx+s1,ry,xc+s2,yc) >= 1) { 
        rx = rx+s1
        xc = xc+s2
    }
    if (check_inner(px,py,rx+s1,ry,xc-s2,yc) >= 1) {
        rx = rx+s1
        xc = xc-s2      
    }
    if (check_inner(px,py,rx,ry+s1,xc,yc+s2) >= 1) {
        ry = ry+s1
        yc = yc+s2
    }
    
    if (check_inner(px,py,rx,ry+s1,xc,yc-s2) >= 1) { 
        ry = ry+s1
        yc = yc-s2      
    }
    if (check_inner(px,py,rx,ry,xc,yc) < 1) {
        rx = rx-2*step-0.001
        ry = ry-2*step-0.001
        if (rx < 0.001)
            rx = 0.001
        if (ry < 0.001)
            ry = 0.001
       ## print ("reduce")
    }

    if (rx == rxprev && ry == ryprev && xcprev == xc && ycprev== yc  && check_inner(px,py,rx,ry,xc,yc) >= 1) { 
        jxy = get_jumperpoint(px,py,step*5,rx,ry,xc,yc)
        jx <- jxy[1]
        jy <- jxy[2]
        ##print ("jump %s %s -> %s %s" % (xc,yc,xc+jx,yc+jy))
        xc = xc+jx
        yc = yc+jy
    }

    if (rx*ry > maxarea) {
        bestiter = iterat
        bestx = xc
        besty = yc
        besta = rx
        bestb = ry
        maxarea = rx*ry
        }

    
}
    return(list(center=c(bestx,besty),radius.x=besta,radius.y=bestb, maxarea=maxarea))
}

#' Inscribe the maximum ellipse into any arbitrary 2D polygon
#'
#' Inscribe the maximum ellipse into any arbitrary 2D polygon
#' @param poly k x 2 matrix containing ordered coordinates forming the polygon. If outline is not closed, the function will close it.
#' @param step stepsize
#' @param iters integer: number of iterations to run
#' 
#' @return 
#' \item{center}{ center of ellipse}
#' \item{radius.x}{ x-dim of ellipse}
#' \item{radius.y}{ y-dim of ellipse}
#' \item{maxarea}{area of ellipse}
#' @examples
#' require(shapes)
#' require(DescTools)
#' poly <- gorf.dat[c(1,6:8,2:5),,1]
#' \dontrun{
#' myellipse <- inscribeEllipse(poly,iters = 200)
#' plot(poly,asp=1)
#' lines(rbind(poly,poly[1,]))
#' DrawEllipse(x=myellipse$center[1],y=myellipse$center[2],radius.x=myellipse$radius.x,
#'             radius.y = myellipse$radius.y,col="red")
#' } 
#' @export
inscribeEllipse <- function(poly,step=0.3,iters=999) {
    if (sum(abs(poly[1,]-poly[nrow(poly),])) != 0)
        poly <- rbind(poly,poly[1,])
    px_old = poly[,1]
    py_old = poly[,2]
    init_point = colMeans(poly)
    init_radius = step
    out <- .Call("inscribeEllipseCpp",poly,step,iters,init_point)
    out$maxarea <- out$maxarea*pi
    return(out)
}


makeRotMat2d <- function(theta) {
    ct <- cos(theta)
    st <- sin(theta)
    out <- matrix(c(ct,st,-st,ct),2,2)
    return(out)
}


#' Inscribe the maximum ellipse into any arbitrary 2D polygon including rotations
#'
#' Inscribe the maximum ellipse into any arbitrary 2D polygon including rotations
#' @param poly k x 2 matrix containing ordered coordinates forming the polygon
#' @param step stepsize
#' @param iters integer: number of iterations to run
#' @param rotsteps integer: number rotational steps 
#' @return 
#' \item{center}{ center of ellipse}
#' \item{radius.x}{ x-dim of ellipse}
#' \item{radius.y}{ y-dim of ellipse}
#' \item{maxarea}{area of ellipse}
#' \item{theta}{angle of optimal rotation}
#' \item{polyRot}{k x 2 matrix of cooridnates rotated around barycenter to maximize ellipse area}
#' @examples
#' require(shapes)
#' require(DescTools)
#' poly <- gorf.dat[c(1,6:8,2:5),,1]
#' \dontrun{
#' myellipse <- inscribeEllipseRot(poly,iters = 999,rotsteps=10)
#' plot(poly,asp=1)
#' lines(rbind(poly,poly[1,]))
#' DrawEllipse(x=myellipse$center[1],y=myellipse$center[2],radius.x=myellipse$radius.x,
#'             radius.y = myellipse$radius.y,col="red")
#' } 
#' @export
inscribeEllipseRot <- function(poly,step=0.3,iters=999,rotsteps=45) {

    thetaList <- seq(0,pi,length.out = rotsteps)[-rotsteps]
    rots <- lapply(thetaList,makeRotMat2d)
    polyS <- scale(poly,scale=FALSE)
    polyRot <- lapply(rots,function(x) x <- polyS%*%x)
    polyRot <- lapply(polyRot,function(x) x <- sweep(x,2,-attributes(polyS)$'scaled:center'))
    init_point = colMeans(poly)
    ## bestarea <- 0
    ## bestfit <- NULL
    ## besttheta <- NULL
    ## for (i in 1:length(polyRot)) {
    ##     tt <- inscribeEllipse(polyRot[[i]],step=step,iters = iters)
    ##     if (tt$maxarea > bestarea) {
    ##         bestarea <- tt$maxarea
    ##         bestfit <- tt
    ##         besti <- i
    ##     }
        
    ## }

    bestfit  <-  .Call("inscribeEllipseRotCpp",polyRot,step,iters,init_point)
    besti <- bestfit$bestiter
    bestfit$theta <- thetaList[besti]
    bestfit$polyRot <- polyRot[[besti]]
    bestfit$maxarea <- bestfit$maxarea*pi
   return(bestfit)
}
