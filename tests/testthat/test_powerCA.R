library(CAs)
test_that("powerCA works", {
  expect_error(powerCA(5,10000,7), "no construction for this setting")
  ## for u3 > 0
  expect_equal(dim(powerCA(3,352,2)), c(49,352))
  expect_equal(dim(powerCA(3,300,2)), c(49,300))
  expect_error(cyc(31,3,k=32), "k can be at most 31")
  ## u4 > 0 for one case only, which alone runs about 3 minutes
  ## kept unchecked
})
