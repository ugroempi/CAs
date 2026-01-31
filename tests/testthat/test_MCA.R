library(CAs)
test_that("MCAt works", {
  ## successful, finds optimum
  D <- MCAt(list(levels=c(6,3,2),frequency=c(2,2,2)), t=3, seed=29120)
  expect_equal(dim(D), c(108,6))
  expect_error(MCAt(c("four","two","six"), 3),
    "if not a list or a data.frame, nlevels must be a vector of numbers of levels")
  expect_error(MCAt(list(levels=c(6,6,3,2),frequency=c(1,1,2,2)), t=3),
               "levels must have unique elements")
  expect_error(MCAt(list(levels=c(6,6,3,2),frequency=c(2,2,2)), t=3),
               "the components 'levels' and 'frequency' must have the same length")
  expect_error(MCAt(list(levels=c(2,3,6),frequency=c(2,2,2)), t=3),
               "must be sorted in decreasing order")
})

test_that("MCA2 works", {
  ## successful, finds optimum
  D <- MCA2(list(levels=c(6,3,2),frequency=c(2,2,2)))
  expect_equal(dim(D), c(36,6))
  expect_error(MCAt(c("four","two","six"), 3),
               "if not a list or a data.frame, nlevels must be a vector of numbers of levels")
  expect_error(MCAt(list(levels=c(6,6,3,2),frequency=c(1,1,2,2)), t=3),
               "levels must have unique elements")
  expect_error(MCAt(list(levels=c(6,6,3,2),frequency=c(2,2,2)), t=3),
               "the components 'levels' and 'frequency' must have the same length")
  expect_error(MCAt(list(levels=c(2,3,6),frequency=c(2,2,2)), t=3),
               "must be sorted in decreasing order")
})

test_that("projBoseMCA works", {
  ## successful, finds optimum
  D <- projBoseMCA(list(levels=c(8,6),frequency=c(4,3)))
  expect_equal(dim(D), c(118,7))
  expect_error(projBoseMCA(c("four","two","six"), 3),
               "if not a list or a data.frame, nlevels must be a vector of numbers of levels")
  expect_error(projBoseMCA(list(levels=c(6,3,2),frequency=c(1,2,2))),
               "projBoseMCA not applicable for more than two different numbers of levels")
  expect_error(projBoseMCA(list(levels=c(8,8),frequency=c(4,3))),
               "levels must have unique elements")
  expect_error(projBoseMCA(list(levels=c(6,3),frequency=c(2,2,2))),
               "the components 'levels' and 'frequency' must have the same length")
  expect_error(projBoseMCA(list(levels=c(3,6),frequency=c(2,2))),
               "must be sorted in decreasing order")
})

test_that("N_upper_MCA works", {
  ## successful, finds optimum
  hilf <- N_upper_MCA(list(levels=c(8,6),frequency=c(4,3)))
  expect_equal(unname(hilf), 64)
  expect_error(N_upper_MCA(c("sechs","drei","zwei")),
               "if not a list or a data.frame, \nnlevels must be a vector of numbers of levels")
  expect_error(N_upper_MCA(list(level=c(8,6),frequency=c(4,3))),
               "nlevels must have components 'levels' and 'frequency'")
  expect_error(N_upper_MCA(list(levels=c(6,3,2),frequency=c(2,2))),
               "the components 'levels' and 'frequency' must have the same length")
})

