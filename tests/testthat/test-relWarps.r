context("relative warps")
test_that("relwarps behaves", {
              data(boneData)
              rW.baselince=readRDS("testdata/rW.rds")
              expect_equal(lapply(relWarps(boneLM),abs),lapply(rW.baseline,abs),tol=1e-7)
})
