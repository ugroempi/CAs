library(testthat)
library(CAs)

## the old (non-cafast) code paths must remain usable
test_that("maxconstant works with nocafast=TRUE", {
  Dpaley <- paleyCA(3,8,2)
  expect_equal(dim(maxconstant(Dpaley, nocafast=TRUE)), c(12,8))
})

test_that("flexpos and markflex work with nocafast=TRUE", {
  D <- rbind(paleyCA(3,8,2), 0)
  fp <- flexpos(D, 2, nocafast=TRUE)
  expect_true(is.matrix(fp))
  expect_equal(dim(fp), dim(D))
  expect_true(is.matrix(suppressMessages(markflex(D, 2, nocafast=TRUE))))
})

test_that("postopNCK works with nocafast=TRUE", {
  D <- rbind(paleyCA(3,8,2), 0)
  aus <- suppressMessages(postopNCK(D, 2, outerRetry=2, innerRetry=1, seed=1, nocafast=TRUE))
  expect_true(nrow(aus) <= nrow(D))
})

test_that("CA_to_PCA works with nocafast=TRUE", {
  D <- paleyCA(3,8,2)
  aus <- suppressMessages(CA_to_PCA(D, nocafast=TRUE))
  expect_s3_class(aus, "ca")
  expect_equal(dim(aus), dim(D))
  expect_false(is.null(attr(aus, "PCAstatus")))
  aus2 <- suppressMessages(CA_to_PCA(D, tryhard=TRUE, nocafast=TRUE))
  expect_s3_class(aus2, "ca")
})
