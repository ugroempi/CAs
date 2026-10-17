library(CAs)

test_that("OD2 uses prime-power, exceptional, and generic constructions", {
  prime_power <- CAs:::OD2(3)
  exceptional6 <- CAs:::OD2(6)
  exceptional10 <- CAs:::OD2(10)
  generic <- CAs:::OD2(8)

  expect_identical(dim(prime_power), c(6L, 3L))
  expect_identical(ncol(exceptional6), 3L)
  expect_identical(ncol(exceptional10), 3L)
  expect_identical(dim(generic), c(56L, 2L))
  expect_true(all(generic >= 0 & generic < 8))
})

test_that("OD3 covers prime-power and full-factorial constructions", {
  prime_power <- CAs:::OD3(3)
  non_prime_power <- CAs:::OD3(6)
  smallest <- CAs:::OD3(1)

  expect_identical(ncol(prime_power), 4L)
  expect_true(all(prime_power %in% 0:3))
  expect_identical(dim(non_prime_power), c(210L, 3L))
  expect_true(all(apply(non_prime_power, 1, function(row) length(unique(row))) == 3L))
  expect_identical(dim(smallest), c(6L, 3L))
  expect_true(all(apply(smallest, 1, function(row) length(unique(row))) == 2L))
})

test_that("cyclic DCA ingredients and their row counts cover each branch", {
  expect_identical(CAs:::circularDCA(7), matrix(7, 1, 1))
  expect_identical(CAs:::circularDCA(1:4),
                   rbind(1:4, c(2:4, 1), c(3:4, 1:2), c(4, 1:3)))

  expect_identical(CAs:::N_DCA(6, 4), 7)
  expect_identical(CAs:::N_DCA(10, 5), 11)
  expect_identical(CAs:::N_DCA(18, 4), 19)
  expect_true(is.na(CAs:::N_DCA(12, 3)))
  expect_equal(CAs:::N_DCA(3, 3), 3)
  expect_equal(CAs:::N_DCA(3, 4), 4)
  expect_equal(CAs:::N_DCA(3, 5), 5)
  expect_equal(CAs:::N_DCA(2, 3), 3)
})

test_that("DCA constructs exceptional, prime-power, and product designs", {
  expect_identical(dim(CAs:::DCA(6, 4)), c(7L, 4L))
  expect_identical(dim(CAs:::DCA(10, 5)), c(11L, 5L))
  expect_identical(dim(CAs:::DCA(18, 4)), c(19L, 4L))
  expect_identical(dim(CAs:::DCA(3, 3)), c(3L, 3L))
  expect_identical(dim(CAs:::DCA(3, 4)), c(4L, 4L))
  expect_identical(dim(CAs:::DCA(3, 5)), c(4L, 5L))
  expect_identical(dim(CAs:::DCA(2, 5)), c(4L, 5L))
  expect_error(CAs:::DCA(12, 3))
})
