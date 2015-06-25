context("relative warps")
test_that("relwarps behaves", {
              data(boneData)
              rW.baseline=readRDS("testdata/rW.rds")
              rwtest <- relWarps(boneLM)
              expect_equal(lapply(rwtest,abs),lapply(rW.baseline,abs),tol=1e-7)
})
