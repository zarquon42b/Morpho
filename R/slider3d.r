slider3d<-function(dat.array,SMvector,outlines=NULL,surp=NULL,sur.path="sur",sur.name=NULL,ignore=NULL,sur.type="ply",clean.init=FALSE,tol=1e-05,deselect=FALSE,inc.check=TRUE,recursive=TRUE,iterations=0,initproc=FALSE,pairedLM=0)

{
  if (iterations == 0)
    {
      iterations <- 1e10
    }
  
  if (is.null(outlines) && is.null(surp))	
    {stop("nothing to slide")}
  
  
  n<-dim(dat.array)[3]
  k<-dim(dat.array)[1]
  m<-dim(dat.array)[2]
  
  if (pairedLM[1]!=0 && is.vector(pairedLM)) # check if there are only 2 symmetric lms
    {pairedLM<-t(as.matrix(pairedLM))
   }
  
### update indexing for after ignored landmarks are removed ###	
  if (!is.null(ignore))
    {li<-length(ignore)
     lm.old<-c(1:k)[-ignore]
     mat.ptr<-matrix(c(1:(k-li),lm.old),k-li,2)
		print(mat.ptr)
                                        #print (mat.ptr,li)
     ptr<-function(xo)	### define pointer function for indexing
       {
         if (length(which(ignore %in% xo))!= 0)
           {xo<-xo[-which(xo %in% ignore)]
          }
         
         for (i in 1:(k-li))
           {xo[which(xo==mat.ptr[i,2])]<-mat.ptr[i,1]
            
          }
         return(xo)
       }
				
     if (!is.null(outlines)) ### update outline indices
       {outlines<-lapply(outlines,ptr)}
     if (!is.null(surp)) 	### update surface indices
       {surp<-ptr(surp)}
     if (!is.null(SMvector)) 
       {SMvector<-ptr(SMvector)### of fixed/sliding definition
      }	
     
     if (pairedLM[1]!=0)	### update paired landmarks indices
       {count<-0
        del<-NULL
        for (i in 1:dim(pairedLM)[1])	
          {if (length(which(ignore %in% pairedLM[i,]))!=0)
             {count<-count+1
              del[count]<-i
            }
         }
        
        pairedLM<-pairedLM[-del,]
        if (is.vector(pairedLM))
          {pairedLM<-t(as.matrix(pairedLM))
         }
        if (dim(pairedLM)[1]==0)
          {pairedLM<-0
         }
        else
          {pairedLM<-apply(pairedLM,2,ptr)
           if (is.vector(pairedLM))
             {pairedLM<-t(as.matrix(pairedLM))
            }
         }
      }
     dat.array<-dat.array[-ignore,,]
     k<-dim(dat.array)[1]
   }
  vn.array<-dat.array
  if(length(sur.name)==0)
    {	
      sur.name<-dimnames(dat.array)[[3]]
      sur.name<-paste(sur.path,"/",sur.name,".",sur.type,sep="")
    }
  if (clean.init==TRUE)
    {clean.mesh(dat.array=dat.array,sur.path=sur.path,sur.name=sur.name,sur.type=sur.type)}
  
  p1<-10^12
  
  ini<-procOPA(dat.array[,,1],dat.array[,,2],reflect=TRUE) # create mean between first tow configs to avoid singular BE Matrix
  mshape<-(ini$Ahat+ini$Bhat)/2
      
  cat(paste("Points will be initially projected onto surfaces","\n","-------------------------------------------","\n"))
  {for (j in 1:n)
     {projBack(dat.array[,,j],sur.name[j])
      a<-read.table("out_cloud.ply",skip=14,sep=" ")
      vs<-as.matrix(a[,1:3])
      vn<-as.matrix(a[,4:6])
      dat.array[,,j]<-vs
      vn.array[,,j]<-vn
      unlink("out_cloud.ply") #clean up
    }
   cat(paste("\n","-------------------------------------------","\n"),"Projection finished","\n","-------------------------------------------","\n")
 }
  
  
  
  if(initproc==TRUE) # perform proc fit before sliding
    {	cat("Inital procrustes fit ...")	
        procini<-procGPA(dat.array,pcaoutput=FALSE,scale=TRUE,reflect=TRUE,distances=FALSE)
        mshape<-procini$mshape
      }
  dataslide<-dat.array
  
  if (pairedLM[1]!=0)# create symmetric mean to get rid of assymetry along outline after first relaxation
    {
      Mir<-diag(c(-1,1,1))
      A<-mshape
      Amir<-mshape%*%Mir
      Amir[c(pairedLM),]<-Amir[c(pairedLM[,2:1]),]
      symproc<-procOPA(A,Amir)
      mshape<-(A+Amir)/2
    }
  cat(paste("Start sliding...","\n","-------------------------------------------","\n"))
  
                                        # calculation for a defined number of iterations
  
  count<-1
  while (p1>tol && count <= iterations)
    { 	dataslide_old<-dataslide
        mshape_old<-mshape
        cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated  
        
        if (recursive==TRUE)      # slided Semilandmarks are used in next iteration step
          { dat.array<-dataslide}
        if (m==3)
          {L<-CreateL(mshape)}
        else 
          {L<-CreateL2D(mshape)} 
        for (j in 1:n)
          { 
            U<-calcTang_U_s(dat.array[,,j],vn.array[,,j],SMvector=SMvector,outlines=outlines,surface=surp,deselect=deselect)
            		dataslide[,,j]<-calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m)$Gamatrix
            
###projection onto surface
            projBack(dataslide[,,j],sur.name[j])
            a<-read.table("out_cloud.ply",skip=14,sep=" ")
            vs<-as.matrix(a[,1:3])
            vn<-as.matrix(a[,4:6])
            dataslide[,,j]<-vs
            vn.array[,,j]<-vn	
            unlink("out_cloud.ply") #clean up
          }
        cat("estimating sample mean shape...")
        proc<-procGPA(dataslide,scale=TRUE,reflect=TRUE,pcaoutput=FALSE,distances=FALSE)
        mshape<-proc$mshape
        if (pairedLM[1]!=0)# create symmetric mean to get rid of assymetry along outline after first relaxation
          {
            Mir<-diag(c(-1,1,1))
            A<-mshape
            Amir<-mshape%*%Mir
            Amir[c(pairedLM),]<-Amir[c(pairedLM[,2:1]),]
            symproc<-procOPA(A,Amir)
      			mshape<-(A+Amir)/2
          }          	
        p1_old<-p1		
        p1<-sum(diag(crossprod((mshape_old/cSize(mshape_old))-(mshape/cSize(mshape)))))
        
### check for increasing convergence criterion ###				
        if (inc.check)
          {
            if (p1 > p1_old)
              {
                dataslide<-dataslide_old
                cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
                p1<-0
              } 
            else
              {
                cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
                count<-count+1         
              }
          }
        else
          {
            cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
            count<-count+1         
          }
      }


  return(list(dataslide=dataslide,vn.array=vn.array))
}
