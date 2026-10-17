library(CAs)

test_that("Galois field arithmetic helpers agree with their lookup tables", {
  gf <- lhs::create_galois_field(4)

  expect_identical(CAs:::int2poly(0:3, gf), gf$poly)
  expect_identical(CAs:::gf_prod(0:3, 3:0, gf),
                   gf$times[cbind(1:4, 4:1)])
  expect_identical(CAs:::gf_sum(0:3, 3:0, gf),
                   gf$plus[cbind(1:4, 4:1)])

  reduced <- CAs:::gf_minus(c(0, 1, 2, 3, 4), c(0, 1, 2, 3, 5), gf)
  expect_identical(CAs:::gf_sum(reduced, c(0, 1, 2, 3, 1), gf),
                   c(0L, 1L, 2L, 3L, 0L))

  summands <- list(c(0, 1, 2), c(1, 2, 3), c(2, 3, 0))
  expect_identical(CAs:::gf_sum_list(summands, gf),
                   CAs:::gf_sum(CAs:::gf_sum(summands[[1]], summands[[2]], gf),
                                summands[[3]], gf))
  factors <- list(c(0, 1, 2), c(3, 2, 1), c(1, 3, 2))
  expect_identical(CAs:::gf_prod_list(factors, gf),
                   CAs:::gf_prod(CAs:::gf_prod(factors[[1]], factors[[2]], gf),
                                 factors[[3]], gf))

  expect_identical(CAs:::gf_pow(c(0, 1, 2), c(0, 1, 2), gf),
                   c(1, 1, CAs:::gf_prod(2, 2, gf)))
  expect_identical(CAs:::gf_pow(c(0, 1, 2), c(3, 2, 1), gf, checks = FALSE),
                   c(CAs:::gf_prod_list(rep(0, 3), gf),
                     CAs:::gf_prod_list(rep(1, 2), gf), 2))
})

test_that("Galois field matrix operations handle checks and determinants", {
  gf <- lhs::create_galois_field(2)
  left <- matrix(c(1, 0, 1, 1), nrow = 2)
  right <- matrix(c(1, 1, 0, 1), nrow = 2)

  expect_identical(CAs:::gf_matmult(left, right, gf),
                   matrix(c(0L, 1L, 1L, 1L), nrow = 2))
  expect_identical(CAs:::gf_matmult(left, right, gf, checks = FALSE),
                   CAs:::gf_matmult(left, right, gf))
  expect_identical(CAs:::gf_det(diag(1, 3), gf), 1L)
  expect_identical(CAs:::gf_det(matrix(1, 2, 2), gf), 0L)
  expect_error(CAs:::gf_matmult(matrix(0, 1, 1), matrix(0, 2, 1), gf))
  expect_error(CAs:::gf_sum_list(list(0:1, 0:0), gf))
  expect_error(CAs:::gf_prod_list(list(0:1, 0:0), gf))
  expect_error(CAs:::gf_pow(1, -1, gf))
  expect_error(CAs:::gf_pow(1, 1.5, gf))
})

test_that("Bush construction and alternate Galois fields are available", {
  expect_identical(dim(CAs:::createBush3PowerOf2(2)), c(8L, 4L))
  expect_identical(dim(CAs:::createBush3PowerOf2(4)), c(64L, 6L))
  expect_error(CAs:::createBush3PowerOf2(3))

  gf8 <- CAs:::mygf(8)
  gf8_2006 <- CAs:::mygf(8, "2006")
  gf8_2009 <- CAs:::mygf(8, "2009")
  expect_equal(gf8$times, gf8_2006$times)
  expect_equal(gf8_2006$times, gf8_2009$times)

  gf9_2006 <- CAs:::mygf(9, "2006")
  gf9_2009 <- CAs:::mygf(9, "2009")
  expect_false(identical(gf9_2006$times, gf9_2009$times))
  expect_identical(dim(gf9_2006$plus), c(9L, 9L))
  expect_identical(length(gf9_2009$inv), 8L)

  for (q in c(16, 25)) {
    gf <- CAs:::mygf(q, "2006")
    expect_identical(dim(gf$times), as.integer(c(q, q)))
    expect_true(all(gf$times[1, ] == 0))
    expect_true(all(gf$times[, 1] == 0))
  }
  expect_error(CAs:::mygf(7))
  expect_error(CAs:::mygf(16, "2009"))
})
