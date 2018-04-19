context("CVA")
test_that("CVA behaves", {
    data(iris)
    vari <- iris[,1:4]
    facto <- iris[,5]
    cva.1 <- CVA(vari, groups=facto)
    CVA.baseline=readRDS("testdata/CVA.rds")
    expect_equal(abs(cva.1$CV),abs(CVA.baseline$CV),tol=0.001,check.attributes=FALSE)
})
