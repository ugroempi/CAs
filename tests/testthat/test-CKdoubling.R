library(CAs)
test_that("CK_doubling works", {
  ## with and without inputs
  ## with and without check
  ## correct and incorrect inputs
  D3correct <- paleyCA(3, 11)
  D3wrong <- KSK(11)
  D3correct1 <- D3correct + 1
  D3wrong1 <- D3wrong + 1
  D2correct <- D3wrong
  D2correct1 <- D3wrong1
  D2wrong <- D2correct[-1,]

  ## start0 both states and different starts (always error)
  ## all specified with start0=TRUE, no check
  D <- CK_doubling(D3correct, D2correct)
  expect_equal(dim(D), c(19, 22))
  expect_equal(max(D), 1)
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## all specified with start0=FALSE, no check
  expect_error(D <- CK_doubling(D3correct1, D2correct1))
  D <- CK_doubling(D3correct1, D2correct1, start0=FALSE)
  expect_equal(dim(D), c(19, 22))
  expect_equal(max(D), 2)
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## all specified with mixed start0, error expected
  expect_error(D <- CK_doubling(D3correct1, D2correct))
  expect_error(D <- CK_doubling(D3correct1, D2correct, start0=FALSE))
  expect_error(D <- CK_doubling(D3correct, D2correct1))
  expect_error(D <- CK_doubling(D3correct, D2correct1, start0=FALSE))

  ## check CK_doublingCA with k=22 and v=2
  D <- CK_doublingCA(22, 2)
  expect_equal(dim(D), c(19, 22))
  expect_equal(max(D), 1)
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## factors D3 correct
  ## D2 correct
  ## check
  ## 8 cases, of which one was already covered
  ## one case with valid result
  D <- CK_doubling(D3correct, D2correct, check=TRUE)
  expect_equal(dim(D), c(19, 22))
  ## six cases with errors or inadequate results
  expect_error(D <- CK_doubling(D3wrong, D2correct, check=TRUE))
  expect_error(D <- CK_doubling(D3correct, D2wrong, check=TRUE))
  D <- CK_doubling(D3wrong, D2correct)
  expect_false(caverify::ca_verify(D, 3)$covered)
  D <- CK_doubling(D3correct, D2wrong)
  expect_false(caverify::ca_verify(D, 3)$covered)
  D <- CK_doubling(D3correct, D2wrong)
  expect_false(caverify::ca_verify(D, 3)$covered)
  D <- CK_doubling(D3wrong, D2wrong)
  expect_false(caverify::ca_verify(D, 3)$covered)
})
