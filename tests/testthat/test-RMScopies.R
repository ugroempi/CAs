library(CAs)
test_that("RMScopies works", {
  ## base construction (Raaphorst-Moura-Stevens 2014)
  D <- RMScopies(3)
  expect_equal(dim(D), c(53, 13))
  expect_true(caverify::ca_verify(D, 3)$covered)
  D <- RMScopies(5)
  expect_equal(dim(D), c(249, 31))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## one shifted copy
  D <- RMScopies(5, 1)
  expect_equal(dim(D), c(324, 32))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## several shifted copies (layers)
  for (m in 2:7){
    D <- RMScopies(5, m)
    expect_equal(nrow(D), N_RMScopiesCA(3, 31 + m, 5))
    expect_equal(ncol(D), 31 + m)
    expect_true(caverify::ca_verify(D, 3)$covered)
  }
  D <- RMScopies(7, 8)
  expect_equal(dim(D), c(980, 65))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## prime powers (Galois field arithmetic from lhs)
  D <- RMScopies(4)
  expect_equal(dim(D), c(127, 21))
  expect_true(caverify::ca_verify(D, 3)$covered)
  D <- RMScopies(4, 1)
  expect_equal(dim(D), c(159, 22))
  expect_true(caverify::ca_verify(D, 3)$covered)
  for (m in c(2, 4, 5, 7)){
    D <- RMScopies(4, m)
    expect_equal(nrow(D), N_RMScopiesCA(3, 21 + m, 4))
    expect_true(caverify::ca_verify(D, 3)$covered)
  }
  D <- RMScopiesCA(74, 8)
  expect_equal(dim(D), c(1407, 74))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## interface function
  E <- RMScopiesCA(37, 5)
  expect_s3_class(E, "ca")
  expect_equal(dim(E), c(350, 37))
  expect_equal(attr(E, "t"), 3)
  expect_equal(dim(RMScopiesCA(20, 5)), c(249, 20))

  ## run sizes
  expect_equal(N_RMScopiesCA(3, 31, 5), 249)
  expect_equal(N_RMScopiesCA(3, 32, 5), 324)
  expect_equal(N_RMScopiesCA(3, 36, 5), 325)
  expect_equal(N_RMScopiesCA(3, 37, 5), 350)
  expect_equal(N_RMScopiesCA(3, 58, 7), 930)
  expect_true(is.na(N_RMScopiesCA(4, 32, 5)))
  expect_true(is.na(N_RMScopiesCA(3, 32, 6)))
  expect_equal(N_RMScopiesCA(3, 81, 8), 1408)
  expect_equal(N_RMScopiesCA(3, 100, 9), 2025)
  expect_true(is.na(N_RMScopiesCA(3, 120, 10)))

  ## number of columns for a given run size
  expect_true(is.na(k_RMScopiesCA(3, 248, 5)))
  expect_equal(k_RMScopiesCA(3, 249, 5), 31)
  expect_equal(k_RMScopiesCA(3, 324, 5), 32)
  expect_equal(k_RMScopiesCA(3, 349, 5), 36)
  expect_equal(k_RMScopiesCA(3, 350, 5), 37)
  expect_equal(k_RMScopiesCA(3, 3751, 11), 144)
  expect_true(is.na(k_RMScopiesCA(4, 350, 5)))
  expect_true(is.na(k_RMScopiesCA(3, 350, 6)))
  ## k_ and N_ work with m=0 outside of data.frame
  ##      and tested_only=TRUE
  expect_equal(k_RMScopiesCA(3, 24333, 23), 553)
  expect_equal(N_RMScopiesCA(3, 400, 23), 24333)
  ## k_ and N_ agree
  for (v in c(7, 8, 9)) for (N in c(930, 1407, 1408, 2025, 3000)){
    k <- k_RMScopiesCA(3, N, v)
    if (is.na(k)) next
    expect_lte(N_RMScopiesCA(3, k, v), N)
    if (k < 2*(v^2 + v + 1)) expect_gt(N_RMScopiesCA(3, k + 1, v), N)
  }

  ## tested parameters
  expect_s3_class(rmsCOPIEScat, "data.frame")
  expect_true(all(c("q", "m", "k", "N", "verifier") %in% names(rmsCOPIEScat)))
  expect_false(any(duplicated(rmsCOPIEScat[, c("q", "m")])))
  for (r in which(rmsCOPIEScat$q <= 5 & rmsCOPIEScat$m %in% c(0, 1, 2, 5)))
    expect_equal(nrow(RMScopies(rmsCOPIEScat$q[r], rmsCOPIEScat$m[r])),
                 rmsCOPIEScat$N[r])
  ## untested combinations: refused unless explicitly requested
  expect_error(RMScopies(11, 50))
  expect_true(is.na(N_RMScopiesCA(3, 183, 11)))
  expect_error(RMScopiesCA(183, 11))
  D <- RMScopies(11, 50, tested_only = FALSE)
  expect_equal(nrow(D), N_RMScopiesCA(3, 183, 11, tested_only = FALSE))
  ## between tested values, the smallest tested m' >= m is used
  expect_equal(N_RMScopiesCA(3, 390, 19), 19855)
  expect_equal(N_RMScopiesCA(3, 140, 11), 3751)

  ## invalid input
  expect_error(RMScopies(6))
  expect_error(RMScopies(12))
  expect_error(RMScopies(2, 1))
  expect_true(is.na(N_RMScopiesCA(3, 8, 2)))
  expect_error(RMScopies(5, 32))
  expect_error(RMScopiesCA(32, 6))
})
