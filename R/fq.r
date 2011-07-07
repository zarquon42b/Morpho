fq <- function(x=NULL)
  {
    image <- FALSE
    history=FALSE
    
    if (is.null(x) || ! x %in% c(0:3))
      {
        ans <- 4
      }
    else
      { ans <- x
      }
    while (! ans%in% c(0:3))
      {
        ans <- readline(" save nothing: 0 \n save workspace: 1\n save workspace and history:2\n cancel: 3\n")
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
    
    
    if (image)
      {
        save.image()
      }
    if (history)
      {  savehistory()
       }
    if (ans != 3)
      {
        pid <- Sys.getpid()
        system(paste("kill -9 ",pid))
      }
  }
    
