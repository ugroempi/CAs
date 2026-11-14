library(CAs)
test_that("JE", {
  expect_error(N_JE(),
               regexp = "At least one of k and t must be specified")
  expect_error(k_JE(),
               regexp = "At least one of N and t must be specified")
  expect_equal(N_JE(4,6), 21)
  expect_true(is.na(N_JE(4,5)))
  expect_true(is.na(N_JE(4,7)))
  expect_equal(k_JE(4,21), 6)
  expect_true(is.na(k_JE(4,20)))

  temp <- JE(4,6)
  expect_equal(dim(temp), c(21, 6))
  expect_equal(colSums(temp), rep(11L, 6))
  expect_s3_class(temp,"ca")
  expect_equal(attr(temp, "origin"), "Johnson and Entringer")
  expect_equal(attr(temp, "comment"), "optimal for k=t+2")
  expect_equal(attr(temp, "t"), 4L)
  expect_equal(prod(unlist(coverage(temp, 4))), 1)

  expect_message(JE(5,8))
  expect_equal(dim(JE(3, 5)), c(10, 5))
  expect_equal(dim(JE(3, 4)), c(10, 4))
  suppressMessages(expect_equal(dim(JE(3, 7)), c(42, 7)))

  ## ZeroSum
  ## check inputs
  expect_error(N_ZeroSum())
  expect_error(N_ZeroSum(3,4))
  expect_error(N_ZeroSum(3,k=NULL,3))
  expect_error(N_ZeroSum(t=NULL,k=4,3))
  expect_error(k_ZeroSum(3,4))
  expect_error(k_ZeroSum(3,N=NULL,3))
  expect_error(k_ZeroSum(t=NULL,N=4,3))

  ## check N and k return values
  expect_equal(N_ZeroSum(4,5,2), 16)
  expect_true(is.na(N_ZeroSum(4,4,2)))
  expect_true(is.na(N_ZeroSum(4,6,2)))
  expect_equal(k_ZeroSum(4,16,2), 5)
  expect_equal(k_ZeroSum(4,20,2), 5)
  expect_true(is.na(k_ZeroSum(4,15,2)))

  ## check arrays
  temp <- ZeroSum(4,5,3)
  expect_equal(dim(temp), c(81, 5))
  expect_equal(colSums(temp), rep(81L, 5))
  expect_s3_class(temp,"ca")
  expect_equal(attr(temp, "origin"), "ZeroSum")
  expect_equal(attr(temp, "comment"), "optimal for k=t+1")
  expect_equal(attr(temp, "t"), 4L)
  expect_equal(prod(unlist(coverage(temp, 4))), 1)
  ## also works in general
  temp <- ZeroSum(4,5,6)
  expect_equal(prod(unlist(coverage(temp, 4))), 1)

  expect_message(ZeroSum(5,8,2))
  suppressMessages(expect_equal(dim(ZeroSum(5, 8, 2)), c(128, 8)))
  expect_equal(dim(ZeroSum(3, 4, 8)), c(512, 4))
})
