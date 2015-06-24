context("relative warps")
test_that("relwarps behaves", {
              data(boneData)
              rW.baseline=readRDS("testdata/rW.rds")
              expect_equal(relWarps(boneLM),rW.baseline)
})
