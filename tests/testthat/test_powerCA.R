library(CAs)
test_that("powerCA works", {
  expect_error(powerCA(5,10000,7), "no construction for this setting")
  ## for u3 > 0
  aus <- powerCA(3,352,2)
  expect_equal(dim(aus), c(49,352))
  expect_true(caverify::ca_verify(aus, 3)$covered)
  ## u4 > 0 for one case only, which alone runs about 3 minutes
  ## kept unchecked
})
