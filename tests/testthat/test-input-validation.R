library(CAs)
library(testthat)

# Tests for input validation error messages

# ========================================
# coverage() function tests
# ========================================
test_that("coverage: validates D parameter", {
  expect_error(
    coverage(D=NULL, t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    coverage(D="not_a_matrix", t=2),
    "D must be a matrix or data.frame"
  )
  expect_error(
    coverage(D=list(1,2,3), t=2),
    "D must be a matrix or data.frame"
  )
})

test_that("coverage: validates t parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)

  expect_error(
    coverage(D),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    coverage(D, t="2"),
    "t must be a single positive integer"
  )
  expect_error(
    coverage(D, t=c(2,3)),
    "t must be a single positive integer"
  )
  expect_error(
    coverage(D, t=0),
    "t must be a single positive integer"
  )
  expect_error(
    coverage(D, t=-1),
    "t must be a single positive integer"
  )
  expect_error(
    coverage(D, t=2.5),
    "t must be a single positive integer"
  )
  expect_error(
    coverage(D, t=5),
    "t \\(5\\) cannot be larger than the number of columns in D \\(2\\)"
  )
})

test_that("coverage: validates isInteger parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)

  expect_error(
    coverage(D, t=2, isInteger="TRUE"),
    "isInteger must be a single logical value"
  )
  expect_error(
    coverage(D, t=2, isInteger=c(TRUE, FALSE)),
    "isInteger must be a single logical value"
  )
  expect_error(
    coverage(D, t=2, isInteger=1),
    "isInteger must be a single logical value"
  )
})

test_that("coverage: validates start0 parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)

  expect_error(
    coverage(D, t=2, start0="TRUE"),
    "start0 must be a single logical value"
  )
  expect_error(
    coverage(D, t=2, start0=c(TRUE, FALSE)),
    "start0 must be a single logical value"
  )
})

test_that("coverage: validates parallel parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)

  expect_error(
    coverage(D, t=2, parallel="1"),
    "parallel must be a single positive integer"
  )
  expect_error(
    coverage(D, t=2, parallel=c(1,2)),
    "parallel must be a single positive integer"
  )
  expect_error(
    coverage(D, t=2, parallel=0),
    "parallel must be a single positive integer"
  )
  expect_error(
    coverage(D, t=2, parallel=-1),
    "parallel must be a single positive integer"
  )
  expect_error(
    coverage(D, t=2, parallel=1.5),
    "parallel must be a single positive integer"
  )
})

test_that("coverage: validates verbose parameter", {
  D <- matrix(c(0,0,1,1,0,1,0,1), ncol=2)

  expect_error(
    coverage(D, t=2, verbose="0"),
    "verbose must be a non-negative integer"
  )
  expect_error(
    coverage(D, t=2, verbose=c(0,1)),
    "verbose must be a non-negative integer"
  )
  expect_error(
    coverage(D, t=2, verbose=-1),
    "verbose must be a non-negative integer"
  )
  expect_error(
    coverage(D, t=2, verbose=1.5),
    "verbose must be a non-negative integer"
  )
})

# ========================================
# MCA2() function tests
# ========================================
test_that("MCA2: validates nlevels parameter", {
  expect_error(
    MCA2(),
    "nlevels must be specified"
  )
  expect_error(
    MCA2(nlevels="not_numeric"),
    "nlevels must be a vector of numbers"
  )
  expect_error(
    MCA2(nlevels=numeric(0)),
    "nlevels must have at least one element"
  )
  expect_error(
    MCA2(nlevels=c(2, 3, NA, 4)),
    "nlevels must have non-NA entries only"
  )
  expect_error(
    MCA2(nlevels=c(2, 3, 1)),
    "all elements of nlevels must be integers >= 2"
  )
  expect_error(
    MCA2(nlevels=c(2, 3, 2.5)),
    "all elements of nlevels must be integers >= 2"
  )
})

test_that("MCA2: validates nlevels data.frame format", {
  expect_error(
    MCA2(nlevels=data.frame(x=c(3,2), y=c(2,3))),
    "nlevels must have components 'levels' and 'frequency'"
  )
  expect_error(
    MCA2(nlevels=data.frame(levels=c("3","2"), frequency=c(2,3))),
    "components 'levels' and 'frequency' in nlevels must be numeric"
  )
  expect_error(
    MCA2(nlevels=data.frame(levels=c(3,2,1), frequency=c(2,3,4))),
    "all 'levels' in nlevels must be integers >= 2"
  )
  expect_error(
    MCA2(nlevels=data.frame(levels=c(3,2), frequency=c(0,3))),
    "all 'frequency' values in nlevels must be positive integers"
  )
  expect_error(
    MCA2(nlevels=data.frame(levels=c(2,3), frequency=c(2,3))),
    "levels in nlevels must be sorted in decreasing order"
  )
  expect_error(
    MCA2(nlevels=list(levels=c(3,2), frequency=c(2,3,4))),
    "components 'levels' and 'frequency' must have the same length"
  )
  expect_error(
    MCA2(nlevels=data.frame(levels=c(3,3,2), frequency=c(2,3,4))),
    "nlevels\\$levels must have unique elements"
  )
})

# ========================================
# MCAt() function tests
# ========================================
test_that("MCAt: validates t parameter", {
  expect_error(
    MCAt(nlevels=c(3,3,2)),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    MCAt(nlevels=c(3,3,2), t="3"),
    "t must be a single positive integer"
  )
  expect_error(
    MCAt(nlevels=c(3,3,2), t=c(3,4)),
    "t must be a single positive integer"
  )
  expect_error(
    MCAt(nlevels=c(3,3,2), t=0),
    "t must be a single positive integer"
  )
  expect_error(
    MCAt(nlevels=c(3,3,2), t=2),
    "t must be greater than 2 \\(for t=2, use MCA2 instead\\)"
  )
})

test_that("MCAt: validates nlevels parameter", {
  expect_error(
    MCAt(t=3),
    "nlevels must be specified"
  )
  expect_error(
    MCAt(nlevels="not_numeric", t=3),
    "nlevels must be a vector of numbers"
  )
  expect_error(
    MCAt(nlevels=numeric(0), t=3),
    "nlevels must have at least one element"
  )
  expect_error(
    MCAt(nlevels=c(2, 3, NA, 4), t=3),
    "nlevels must have non-NA entries only"
  )
  expect_error(
    MCAt(nlevels=c(2, 3, 1), t=3),
    "all elements of nlevels must be integers >= 2"
  )
})

test_that("MCAt: validates nlevels data.frame format", {
  expect_error(
    MCAt(nlevels=data.frame(x=c(3,2), y=c(2,3)), t=3),
    "nlevels must have components 'levels' and 'frequency'"
  )
  expect_error(
    MCAt(nlevels=data.frame(levels=c("3","2"), frequency=c(2,3)), t=3),
    "components 'levels' and 'frequency' in nlevels must be numeric"
  )
  expect_error(
    MCAt(nlevels=data.frame(levels=c(3,2,1), frequency=c(2,3,4)), t=3),
    "all 'levels' in nlevels must be integers >= 2"
  )
  expect_error(
    MCAt(nlevels=data.frame(levels=c(3,2), frequency=c(0,3)), t=3),
    "all 'frequency' values in nlevels must be positive integers"
  )
  expect_error(
    MCAt(nlevels=data.frame(levels=c(2,3), frequency=c(2,3)), t=3),
    "levels in nlevels must be sorted in decreasing order"
  )
  expect_error(
    MCAt(nlevels=list(levels=c(3,2), frequency=c(2,3,4)), t=3),
    "components 'levels' and 'frequency' must have the same length"
  )
  expect_error(
    MCAt(nlevels=data.frame(levels=c(3,3,2), frequency=c(2,3,4)), t=3),
    "nlevels\\$levels must have unique elements"
  )
})

# ========================================
# powerCA() function tests
# ========================================
test_that("powerCA: validates t parameter", {
  expect_error(
    powerCA(k=5, v=2),
    "t \\(interaction strength\\) must be specified"
  )
  expect_error(
    powerCA(t="3", k=5, v=2),
    "t must be a single positive integer"
  )
  expect_error(
    powerCA(t=c(3,4), k=5, v=2),
    "t must be a single positive integer"
  )
  expect_error(
    powerCA(t=0, k=5, v=2),
    "t must be a single positive integer"
  )
  expect_error(
    powerCA(t=2.5, k=5, v=2),
    "t must be a single positive integer"
  )
})

test_that("powerCA: validates k parameter", {
  expect_error(
    powerCA(t=3, v=2),
    "k \\(number of columns\\) must be specified"
  )
  expect_error(
    powerCA(t=3, k="5", v=2),
    "k must be a single positive integer"
  )
  expect_error(
    powerCA(t=3, k=c(5,6), v=2),
    "k must be a single positive integer"
  )
  expect_error(
    powerCA(t=3, k=0, v=2),
    "k must be a single positive integer"
  )
})

test_that("powerCA: validates v parameter", {
  expect_error(
    powerCA(t=3, k=5),
    "v \\(number of levels\\) must be specified"
  )
  expect_error(
    powerCA(t=3, k=5, v="2"),
    "v must be a single integer >= 2"
  )
  expect_error(
    powerCA(t=3, k=5, v=c(2,3)),
    "v must be a single integer >= 2"
  )
  expect_error(
    powerCA(t=3, k=5, v=1),
    "v must be a single integer >= 2"
  )
})

test_that("powerCA: validates type parameter", {
  expect_error(
    powerCA(t=3, k=5, v=2, type=123),
    "type must be a single character string"
  )
  expect_error(
    powerCA(t=3, k=5, v=2, type=c("CT","other")),
    "type must be a single character string"
  )
  expect_error(
    powerCA(t=3, k=5, v=2, type="XYZ"),
    "currently only type='CT' is implemented"
  )
})

test_that("powerCA: validates construction availability", {
  expect_error(
    powerCA(t=10, k=5, v=2),
    "no construction for this setting \\(t=10, k=5, v=2\\)"
  )
})

# ========================================
# KSK() function tests
# ========================================
test_that("KSK: validates k and N parameters", {
  expect_error(
    KSK(),
    "At least one of k and N must be specified"
  )
  expect_error(
    KSK(k="10"),
    "k must be a single positive integer"
  )
  expect_error(
    KSK(k=c(10,20)),
    "k must be a single positive integer"
  )
  expect_error(
    KSK(k=0),
    "k must be a single positive integer"
  )
  expect_error(
    KSK(k=-5),
    "k must be a single positive integer"
  )
  expect_error(
    KSK(k=10.5),
    "k must be a single positive integer"
  )
  expect_error(
    KSK(N="10"),
    "N must be a single integer >= 2"
  )
  expect_error(
    KSK(N=c(10,20)),
    "N must be a single integer >= 2"
  )
  expect_error(
    KSK(N=1),
    "N must be a single integer >= 2"
  )
  expect_error(
    KSK(N=-5),
    "N must be a single integer >= 2"
  )
  expect_error(
    KSK(N=10.5),
    "N must be a single integer >= 2"
  )
})

# ========================================
# CAEX() function tests
# ========================================
test_that("CAEX: validates t parameter", {
  expect_error(
    CAEX(k=5, t="2"),
    "t must be a single integer"
  )
  expect_error(
    CAEX(k=5, t=c(2,3)),
    "t must be a single integer"
  )
  expect_error(
    CAEX(k=5, t=3),
    "CAEX requires t=2"
  )
})

test_that("CAEX: validates v parameter", {
  expect_error(
    CAEX(k=5, v="3"),
    "v must be a single integer"
  )
  expect_error(
    CAEX(k=5, v=c(3,4)),
    "v must be a single integer"
  )
  expect_error(
    CAEX(k=5, v=4),
    "CAEX requires v=3"
  )
})

test_that("CAEX: validates k and N parameters", {
  expect_error(
    CAEX(),
    "at least one of k and N must be specified"
  )
  expect_error(
    CAEX(k="5"),
    "k must be a single positive integer"
  )
  expect_error(
    CAEX(k=c(5,6)),
    "k must be a single positive integer"
  )
  expect_error(
    CAEX(k=0),
    "k must be a single positive integer"
  )
  expect_error(
    CAEX(k=-5),
    "k must be a single positive integer"
  )
  expect_error(
    CAEX(N="10"),
    "N must be a single integer"
  )
  expect_error(
    CAEX(N=c(10,11)),
    "N must be a single integer"
  )
  expect_error(
    CAEX(N=5),
    "N must be at least 9"
  )
})

test_that("CAEX: validates maxk1 parameter", {
  expect_error(
    CAEX(k=5, maxk1="TRUE"),
    "maxk1 must be a single logical value"
  )
  expect_error(
    CAEX(k=5, maxk1=c(TRUE, FALSE)),
    "maxk1 must be a single logical value"
  )
  expect_error(
    CAEX(k=5, maxk1=1),
    "maxk1 must be a single logical value"
  )
})

test_that("CAEX: validates k accommodability", {
  expect_error(
    CAEX(k=999999),
    "The requested k=999999 cannot be accommodated"
  )
})

# ========================================
# CA_to_MCA() function tests
# ========================================
test_that("CA_to_MCA: validates D parameter", {
  expect_error(
    CA_to_MCA(cs=1, tolevs=2, t=2),
    "D \\(covering array\\) must be specified"
  )
  expect_error(
    CA_to_MCA(D=data.frame(a=1:3), cs=1, tolevs=2, t=2),
    "D must be a matrix"
  )
  expect_error(
    CA_to_MCA(D=matrix(c("a","b","c","d"), ncol=2), cs=1, tolevs=2, t=2),
    "D must be numeric"
  )
})

test_that("CA_to_MCA: validates cs and tolevs parameters", {
  D <- matrix(c(0,0,0,1,1,1,2,2,2,0,1,2,0,1,2,0,1,2), ncol=3)

  expect_error(
    CA_to_MCA(D=D, tolevs=2, t=2),
    "cs \\(column numbers to modify\\) must be specified"
  )
  expect_error(
    CA_to_MCA(D=D, cs=1, t=2),
    "tolevs \\(target number of levels\\) must be specified"
  )
  expect_error(
    CA_to_MCA(D=D, cs="1", tolevs=2, t=2),
    "cs must be a numeric vector of column numbers"
  )
  expect_error(
    CA_to_MCA(D=D, cs=c(1,0), tolevs=c(2,2), t=2),
    "all elements in cs must be positive integers"
  )
  expect_error(
    CA_to_MCA(D=D, cs=c(1,5), tolevs=c(2,2), t=2),
    "all elements in cs must be valid column numbers"
  )
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs="2", t=2),
    "tolevs must be a numeric vector"
  )
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=1, t=2),
    "all elements in tolevs must be integers >= 2"
  )
  expect_error(
    CA_to_MCA(D=D, cs=c(1,2), tolevs=2, t=2),
    "cs and tolevs must have the same length"
  )
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=5, t=2),
    "tolevs cannot be larger than the current number of levels"
  )
})

test_that("CA_to_MCA: validates t parameter", {
  D <- matrix(c(0,0,0,1,1,1,2,2,2,0,1,2,0,1,2,0,1,2), ncol=3)
  attr(D, "t") <- NULL

  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=2),
    "D does not have a valid attribute 't'"
  )

  attr(D, "t") <- "2"
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=2),
    "D does not have a valid attribute 't'"
  )

  D <- matrix(c(0,0,0,1,1,1,2,2,2,0,1,2,0,1,2,0,1,2), ncol=3)
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=2, t="2"),
    "t must be a single positive integer"
  )
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=2, t=c(2,3)),
    "t must be a single positive integer"
  )
  expect_error(
    CA_to_MCA(D=D, cs=1, tolevs=2, t=0),
    "t must be a single positive integer"
  )
})
