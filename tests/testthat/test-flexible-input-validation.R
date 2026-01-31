library(CAs)
library(testthat)

# Tests for input validation error messages in flexible.R functions

# ========================================
# markflex() function tests
# ========================================
test_that("markflex: validates D parameter", {
  expect_error(
    markflex(D=NULL, t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    markflex(D="not_a_matrix", t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    markflex(D=list(1,2,3), t=2),
    "D must be a matrix or data.frame"
  )
})

test_that("markflex: validates t parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    markflex(D),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    markflex(D, t="2"),
    "t must be a single positive integer"
  )
  expect_error(
    markflex(D, t=c(2,3)),
    "t must be a single positive integer"
  )
  expect_error(
    markflex(D, t=0),
    "t must be a single positive integer"
  )
  expect_error(
    markflex(D, t=-1),
    "t must be a single positive integer"
  )
  expect_error(
    markflex(D, t=2.5),
    "t must be a single positive integer"
  )
  expect_error(
    markflex(D, t=1),
    "t must be at least 2"
  )
  expect_error(
    markflex(D, t=5),
    "t \\(5\\) cannot be larger than the number of columns in D \\(2\\)"
  )
})

test_that("markflex: validates fixrows parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    markflex(D, t=2, fixrows="0"),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    markflex(D, t=2, fixrows=c(0,1)),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    markflex(D, t=2, fixrows=-1),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    markflex(D, t=2, fixrows=1.5),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    markflex(D, t=2, fixrows=10),
    "fixrows \\(10\\) cannot be larger than the number of rows in D \\(4\\)"
  )
})

test_that("markflex: validates verbose parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    markflex(D, t=2, verbose="0"),
    "verbose must be a non-negative integer"
  )
  expect_error(
    markflex(D, t=2, verbose=c(0,1)),
    "verbose must be a non-negative integer"
  )
  expect_error(
    markflex(D, t=2, verbose=-1),
    "verbose must be a non-negative integer"
  )
  expect_error(
    markflex(D, t=2, verbose=1.5),
    "verbose must be a non-negative integer"
  )
})

# ========================================
# flexpos() function tests
# ========================================
test_that("flexpos: validates D parameter", {
  expect_error(
    flexpos(D=NULL, t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    flexpos(D="not_a_matrix", t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    flexpos(D=list(1,2,3), t=2),
    "D must be a matrix or data.frame"
  )
})

test_that("flexpos: validates t parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    flexpos(D),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    flexpos(D, t="2"),
    "t must be a single positive integer"
  )
  expect_error(
    flexpos(D, t=c(2,3)),
    "t must be a single positive integer"
  )
  expect_error(
    flexpos(D, t=0),
    "t must be a single positive integer"
  )
  expect_error(
    flexpos(D, t=-1),
    "t must be a single positive integer"
  )
  expect_error(
    flexpos(D, t=2.5),
    "t must be a single positive integer"
  )
  expect_error(
    flexpos(D, t=1),
    "t must be at least 2"
  )
  expect_error(
    flexpos(D, t=5),
    "t \\(5\\) cannot be larger than the number of columns in D \\(2\\)"
  )
})

# ========================================
# flexprofile() function tests
# ========================================
test_that("flexprofile: validates D parameter", {
  expect_error(
    flexprofile(D=NULL),
    "D must be a matrix or data.frame"
  )
  expect_error(
    flexprofile(D="not_a_matrix"),
    "D must be a matrix or data.frame"
  )
  expect_error(
    flexprofile(D=list(1,2,3)),
    "D must be a matrix or data.frame"
  )
})

# ========================================
# postopNCK() function tests
# ========================================
test_that("postopNCK: validates D parameter", {
  expect_error(
    postopNCK(D=NULL, t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    postopNCK(D="not_a_matrix", t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    postopNCK(D=list(1,2,3), t=2),
    "D must be a matrix or data.frame"
  )
})

test_that("postopNCK: validates t parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    postopNCK(D, t="2"),
    "t must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=c(2,3)),
    "t must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=0),
    "t must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=-1),
    "t must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2.5),
    "t must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=1),
    "t must be at least 2"
  )
  expect_error(
    postopNCK(D, t=5),
    "t \\(5\\) cannot be larger than the number of columns in D \\(2\\)"
  )
})

test_that("postopNCK: validates fixrows parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, fixrows="0"),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, fixrows=c(0,1)),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, fixrows=-1),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, fixrows=1.5),
    "fixrows must be a single non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, fixrows=10),
    "fixrows \\(10\\) cannot be larger than the number of rows in D \\(4\\)"
  )
})

test_that("postopNCK: validates verbose parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, verbose="0"),
    "verbose must be a non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, verbose=c(0,1)),
    "verbose must be a non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, verbose=-1),
    "verbose must be a non-negative integer"
  )
  expect_error(
    postopNCK(D, t=2, verbose=1.5),
    "verbose must be a non-negative integer"
  )
})

test_that("postopNCK: validates outerRetry parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, outerRetry="50"),
    "outerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerRetry=c(50,60)),
    "outerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerRetry=0),
    "outerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerRetry=-1),
    "outerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerRetry=1.5),
    "outerRetry must be a single positive integer"
  )
})

test_that("postopNCK: validates outerMaxnochange parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, outerMaxnochange="10"),
    "outerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerMaxnochange=c(10,20)),
    "outerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerMaxnochange=0),
    "outerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerMaxnochange=-1),
    "outerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, outerMaxnochange=1.5),
    "outerMaxnochange must be a single positive integer"
  )
})

test_that("postopNCK: validates innerRetry parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, innerRetry="10"),
    "innerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerRetry=c(10,20)),
    "innerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerRetry=0),
    "innerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerRetry=-1),
    "innerRetry must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerRetry=1.5),
    "innerRetry must be a single positive integer"
  )
})

test_that("postopNCK: validates innerMaxnochange parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, innerMaxnochange="25"),
    "innerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerMaxnochange=c(25,30)),
    "innerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerMaxnochange=0),
    "innerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerMaxnochange=-1),
    "innerMaxnochange must be a single positive integer"
  )
  expect_error(
    postopNCK(D, t=2, innerMaxnochange=1.5),
    "innerMaxnochange must be a single positive integer"
  )
})

test_that("postopNCK: validates seed parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    postopNCK(D, t=2, seed="123"),
    "seed must be NULL or a single integer"
  )
  expect_error(
    postopNCK(D, t=2, seed=c(123,456)),
    "seed must be NULL or a single integer"
  )
  expect_error(
    postopNCK(D, t=2, seed=1.5),
    "seed must be NULL or a single integer"
  )
})

# ========================================
# uniquecount() function tests
# ========================================
test_that("uniquecount: validates D parameter", {
  expect_error(
    uniquecount(D=NULL, t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    uniquecount(D="not_a_matrix", t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    uniquecount(D=list(1,2,3), t=2),
    "D must be a matrix or data.frame"
  )
})

test_that("uniquecount: validates t parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)
  
  expect_error(
    uniquecount(D),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    uniquecount(D, t="2"),
    "t must be a single positive integer"
  )
  expect_error(
    uniquecount(D, t=c(2,3)),
    "t must be a single positive integer"
  )
  expect_error(
    uniquecount(D, t=0),
    "t must be a single positive integer"
  )
  expect_error(
    uniquecount(D, t=-1),
    "t must be a single positive integer"
  )
  expect_error(
    uniquecount(D, t=2.5),
    "t must be a single positive integer"
  )
  expect_error(
    uniquecount(D, t=1),
    "t must be at least 2"
  )
  expect_error(
    uniquecount(D, t=5),
    "t \\(5\\) cannot be larger than the number of columns in D \\(2\\)"
  )
})
