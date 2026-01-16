library(CAs)
test_that("miscCA", {
  expect_error(miscCA(),
               regexp = '"t"')
  expect_error(miscCA(3),
               regexp = '"k"')
  expect_error(miscCA(3,12),
               regexp = '"v"')
  temp <- miscCA(3,8,2)
  expect_equal(dim(temp), c(12,8))
  expect_equal(prod(unlist(coverage(temp, 3))), 1)

  ## Cohen
  expect_equal(dim(miscCA(2,10,3)), c(14,10))
  ## DoE.base
  expect_equal(dim(temp <- maxconstant(miscCA(2,3,6))), c(36,3))
  expect_equal(length(unique(temp[6,])),1)
  ## Ji-Yin
  expect_equal(dim(temp <- miscCA(3,6,12)), c(1728, 6))
  expect_equal(attr(temp, "PCAstatus")$k1, 4)
  expect_equal(dim(temp <- miscCA(3,6,15)), c(3375, 6))
  expect_equal(attr(temp, "PCAstatus")$k1, 4)
  ## Kokkaly, only for constant rows
  expect_equal(nrow(temp <- miscCA(2,7,5)),30)
  expect_true(is.null(attr(temp, "PCAstatus")))
  expect_true(length(unique(temp[3,]))==1)
  expect_s3_class(temp,"ca")
})
