relaxLM <- function(lm,reference,SMvector,outlines=NULL,surp=NULL,sur.name=NULL,tol=1e-05,deselect=FALSE,inc.check=TRUE,iterations=0)
  {

    k <- dim(lm)[1]
    m <- dim(lm)[2]
    p1<-10^12
    L<-CreateL(reference)
    
    if (iterations == 0)
      {
        iterations <- 1e10
      }


    cat(paste("Points will be initially projected onto surfaces","\n","-------------------------------------------","\n"))
    projBack(lm,sur.name)
    a<-read.table("out_cloud.ply",skip=14,sep=" ")
    vs<-as.matrix(a[,1:3])
    vn<-as.matrix(a[,4:6])
    unlink("out_cloud.ply") #clean up

    count<-1
    while (p1>tol && count <= iterations)
      {
        lm_old <- vs
        cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated
        U<-calcTang_U_s(vs,vn,SMvector=SMvector,outlines=outlines,surface=surp,deselect=deselect)
        dataslido<-calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m)$Gamatrix
        projBack(dataslido,sur.name)
        a<-read.table("out_cloud.ply",skip=14,sep=" ")
        vs<-as.matrix(a[,1:3])
        vn<-as.matrix(a[,4:6])
        unlink("out_cloud.ply") #clean up

        p1_old<-p1
        testproc<-rotonto(lm_old,vs)			   	
        p1<-sum(diag(crossprod((testproc$X/cSize(testproc$X))-(testproc$Y/cSize(testproc$Y)))))### check for increasing convergence criterion ###		
        if (inc.check)
          {
            if (p1 > p1_old)
              {dataslide<-dataslide_old
               cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
               p1<-0
             } 
            else
              {cat(paste("squared distance between iterations:",p1,sep=" "),"\n","-------------------------------------------","\n")
               count<-count+1         
             }
            
          }
      }
    return(vs)
  }
    


