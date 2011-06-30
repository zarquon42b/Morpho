fq <- function(ask=TRUE,image=TRUE,history=FALSE)
  {
    if (ask)
      ans <- 3
      { while (! ans%in% c(0:2))
        {
          ans <- readline(" save nothing: 0 \n save workspace: 1\n save workspace and history:2\n")
        }
        if (ans == 0)
          {image <- FALSE
           history <- FALSE
         }
        else if (ans == 1)
          { image <- TRUE
           history <- FALSE
          }
         else if (ans == 2)
           { image <- TRUE
           history <- TRUE
          }
      }
        
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
