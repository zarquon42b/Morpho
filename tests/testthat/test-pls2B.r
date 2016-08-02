context("2Block PLS")
test_that("pls2B behaves",{
    library(shapes)
    pls.baseline <- readRDS("testdata/pls-baseline.rds")
    proc <- procSym(gorf.dat)
    set.seed(42)
    pls1 <- pls2B(proc$rotated[1:4,,],proc$rotated[5:8,,],same.config=TRUE,rounds=0,mc.cores=1)
    expect_equal(lapply(pls1$svd,abs),lapply(pls.baseline$svd,abs),tol=0.01)
})
