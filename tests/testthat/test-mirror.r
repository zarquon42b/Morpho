context("mirror functions")

test_that("mirror.matrix behaves", {
  boneMir.baseline=structure(c(125.048055055629, 125.337650940873, 121.294041363046, 
    128.549755916167, 117.262978257761, 128.708329677512, 114.015199391521, 
    133.141679085682, 122.716932933575, 123.488377378233, 44.6403187557201, 
    44.9064753775645, 48.309875510538, 48.5435509841916, 55.9361737827255, 
    56.3929693160991, 56.8352669034508, 57.657380761555, 54.4230460742932, 
    40.7425425338623, 76.4985127131541, 80.057848539667, 78.1445857853691, 
    77.6344304495789, 28.8906474576085, 28.8284524490006, 31.1040380007055, 
    31.423573621832, 22.511397330662, 44.668913652422), .Dim = c(10L, 
    3L))
  data(boneData)
  expect_equal(mirror(boneLM[,,1],icpiter=50), boneMir.baseline, tol=1e-6)
  
  skullMir.baseline=readRDS("testdata/skullMir.rds")
  # subsample involves non-deterministic k-means
  set.seed(42)
  expect_equal(mirror(skull_0144_ch_fe.mesh,icpiter=1,subsample = 30,pcAlign=T),
              skullMir.baseline,tol=0.5)
})
