context("sliding semilandmarks")
test_that("slider behaves", {
              data(boneData)
              slider.baseline=readRDS("testdata/slider3d.rds")
              data(nose)
              longnose.mesh <- tps3d(shortnose.mesh,shortnose.lm,longnose.lm,threads=1)
              meshlist <- list(shortnose.mesh,longnose.mesh)
              data <- bindArr(shortnose.lm, longnose.lm, along=3)
              dimnames(data)[[3]] <- c("shortnose", "longnose")
              fix <- c(1:5,20:21)
              surp <- c(1:nrow(shortnose.lm))[-fix]
              slide <- slider3d(data, SMvector=fix, deselect=TRUE,meshlist=meshlist,surp=surp,iterations=1,mc.cores=1,fixRepro=FALSE)$dataslide
              expect_equal(slide,slider.baseline,tol=0.01)
})
