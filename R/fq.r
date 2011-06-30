fq <- function(image=TRUE,history=FALSE)
  {
    if (image)
      {
      save.image()
    }
    if (history)
      {  savehistory()
       }
    pid <- Sys.getpid()
    system(paste("kill -9 ",pid))
  }
