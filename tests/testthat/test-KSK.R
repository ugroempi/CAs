test_that("KSK", {
  expect_error(KSK(),
               regexp = "At least one of k and N must be specified")
  expect_error(KSK(100,4))

  temp <- KSK(100)
  expect_equal(dim(temp), c(10, 100))
  expect_equal(colSums(temp), rep(5L, 100))
  expect_s3_class(temp,"ca")
  expect_equal(attr(temp, "origin"), "Kleitman & Spencer, and Katona")
  expect_equal(attr(temp, "t"), 2L)
  expect_equal(prod(unlist(coverage(temp, 2))), 1)

  expect_warning(expect_error(KSK(N=100)))
  expect_equal(dim(KSK(N=10)), c(10, 126))
  expect_equal(dim(KSK(k=5, N=20)), c(20, 5))

  expect_equal(N_KSK(20), 8)
  expect_equal(k_KSK(12), 462)
})
