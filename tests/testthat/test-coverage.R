library(CAs)
test_that("coverage", {
         perfect <- coverage(bestCA(3,8,4),3)

    ## a simple mixed-level CA of strength 2
    ## that has also a few flexible values
    D <- MCA2(nlevels=c(6,4,3,2,2,2,2,2,2), seed=2481) ## 24 runs
    cD1 <- coverage(D, 2, verbose=1)
    cD2 <- coverage(D, 2, verbose=2)
    cD13 <- coverage(D, 3, verbose=1)
    cD23 <- coverage(D, 3, verbose=2)

    expect_is(cD1, "coverage")
    expect_is(cD2, "coverage")
    expect_is(cD13, "coverage")
    expect_is(cD23, "coverage")

    expect_identical(length(cD1), 7L)
    expect_identical(length(cD2), 10L)
    expect_identical(length(cD13), 7L)
    expect_identical(length(cD23), 10L)

    ##
    expect_identical(cD1[1:4], unclass(perfect))
    expect_identical(cD2[1:4], unclass(perfect))
    expect_identical(cD13$total*332, 179)
    expect_identical(cD13$ave*1728, 1103)
    expect_identical(cD13$min*24, 7)
    expect_identical(cD13$simple*21, 1)

    ##
    expect_true(length(cD1$ncovereds)==36L)
    expect_true(length(cD13$ncovereds)==84L)
    expect_identical(unname(lengths(cD23)), c(rep(1L, 4), rep(84L, 3), 3L*84L, 84L, 84L))
    expect_is(cD23$tabs, "list")

    D <- as.data.frame(D)
    cD1 <- coverage(D, 2, verbose=1)
    cD2 <- coverage(D, 2, verbose=2)
    cD13 <- coverage(D, 3, verbose=1)
    cD23 <- coverage(D, 3, verbose=2)

    expect_is(cD1, "coverage")
    expect_is(cD2, "coverage")
    expect_is(cD13, "coverage")
    expect_is(cD23, "coverage")

    expect_identical(length(cD1), 7L)
    expect_identical(length(cD2), 10L)
    expect_identical(length(cD13), 7L)
    expect_identical(length(cD23), 10L)

    ##
    expect_identical(cD1[1:4], unclass(perfect))
    expect_identical(cD2[1:4], unclass(perfect))
    expect_identical(cD13$total*332, 179)
    expect_identical(cD13$ave*1728, 1103)
    expect_identical(cD13$min*24, 7)
    expect_identical(cD13$simple*21, 1)

    ##
    expect_true(length(cD1$ncovereds)==36L)
    expect_true(length(cD13$ncovereds)==84L)
    expect_identical(unname(lengths(cD23)), c(rep(1L, 4), rep(84L, 3), 3L*84L, 84L, 84L))
    expect_is(cD23$tabs, "list")

#######################################################################
## create test case designs from a
## strength 3 CA with 24 runs:
#######################################################################
##    factor A: t_D=2, v=4; t_D=3, v=3; t_D=4, v=2; t_D=5, v=2
##             list position + 1 is t_D
##    factor B: integer-valued design    ## _noninteger indicates non-integers
##    factor C: start0                   ## _start1 indicates start1 for integers
##    factor D: data frame               ## Ddflist indicates data frame
##    factor E: isInteger of the call (0=FALSE, 1=TRUE)
##    factor F: start0 of the call (0=FALSE, 1=TRUE)
##    factor G: t=t_D or t=t_D+1
##    factor H: parallel=1 or parallel=4
#######################################################################

t_D <- 2:5  ## factor A: t_D,
    ## realized via its element in the respective D list
    ## directly aligns with v
v_D <- c(4,3,2,2)
k <- 7 ## arbitrary, assuming that k is not relevant,
       ## apart from gracefully failing for edge cases
Dlist <- mapply(bestCA, t_D, rep(k,4), v_D)
## entries of Dlist are integer matrices starting with 0,
##    i.e. B TRUE, C TRUE, D FALSE
Dlist_start1 <- lapply(Dlist, function(obj) obj+1)
##         B TRUE, C FALSE, D FALSE
Dlist_noninteger <- lapply(Dlist_start1, function(obj){
  hilf <- obj
  hilf <- letters[hilf]
  dim(hilf) <- dim(obj)
  hilf
})
##         B FALSE, C irrelevant, D FALSE
##                  C may react different in code, i.e., needs to be there
Ddflist <- lapply(Dlist, as.data.frame)
##    i.e. B TRUE, C TRUE, D TRUE
Ddflist_start1 <- lapply(Dlist_start1, as.data.frame)
##         B TRUE, C FALSE, D TRUE
Ddflist_noninteger <- lapply(Dlist_noninteger, as.data.frame)
##         B FALSE, C irrelevant, D FALSE
##                  C may react different in code, i.e., needs to be there
############################################################################

perfect <- coverage(bestCA(3,8,4),3)

    testthat::test_that("CA-based test for coverage inputs (uniform inputs only)", {
    testthat::expect_error(
      coverage(Ddflist_start1[[4]], t=6, isInteger=TRUE, start0=TRUE, parallel=4),
      regexp="start0 is TRUE, but columns do not start at 0"
    )
    testthat::expect_identical(
      coverage(Dlist[[2]], t=4, isInteger=FALSE, start0=FALSE, parallel=4)$simple,
      0)
    testthat::expect_identical(
      coverage(Dlist_start1[[1]], t=2, isInteger=TRUE, start0=FALSE, parallel=4),
      perfect)
    testthat::expect_identical(
      coverage(Ddflist_noninteger[[4]], t=6, isInteger=FALSE, start0=FALSE, parallel=4)$simple,
      0)
    testthat::expect_error(
      coverage(Dlist_noninteger[[3]], t=5, isInteger=TRUE, start0=TRUE, parallel=4),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    testthat::expect_identical(
      coverage(Dlist_noninteger[[2]], t=4, isInteger=FALSE, start0=TRUE, parallel=1)$simple,
      0)
    testthat::expect_identical(
      coverage(Dlist_noninteger[[4]], t=5, isInteger=FALSE, start0=TRUE, parallel=4),
      perfect)
    testthat::expect_identical(
      coverage(Ddflist[[2]], t=3, isInteger=TRUE, start0=TRUE, parallel=4),
      perfect)
## modified isInteger setting for the following row,
## in order to have also error with start0=FALSE
    testthat::expect_error(
      coverage(Ddflist[[3]], t=4, isInteger=TRUE, start0=FALSE, parallel=4),
      regexp="start0 is FALSE, but columns do not start at 1"
    )
    testthat::expect_error(
      coverage(Dlist_start1[[3]], t=4, isInteger=TRUE, start0=TRUE, parallel=1),
      regexp="start0 is TRUE, but columns do not start at 0"
    )
    testthat::expect_error(
      coverage(Ddflist_noninteger[[1]], t=2, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    testthat::expect_identical(
      coverage(Dlist[[1]], t=3, isInteger=FALSE, start0=FALSE, parallel=1)$simple,
      0)
    testthat::expect_identical(
      coverage(Ddflist[[1]], t=3, isInteger=TRUE, start0=TRUE, parallel=1)$simple,
      0)
    testthat::expect_identical(
      coverage(Dlist_noninteger[[3]], t=4, isInteger=FALSE, start0=FALSE, parallel=1),
      perfect)
    testthat::expect_identical(
      coverage(Ddflist_noninteger[[1]], t=2, isInteger=FALSE, start0=TRUE, parallel=4),
      perfect)
    testthat::expect_error(
      coverage(Ddflist_noninteger[[2]], t=4, isInteger=TRUE, start0=TRUE, parallel=4),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    testthat::expect_identical(
      coverage(Ddflist_start1[[2]], t=3, isInteger=FALSE, start0=FALSE, parallel=1),
      perfect)
    testthat::expect_error(
      coverage(Dlist_noninteger[[2]], t=3, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    testthat::expect_error(
      coverage(Ddflist_noninteger[[3]], t=5, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    testthat::expect_identical(
      coverage(Dlist_noninteger[[1]], t=3, isInteger=FALSE, start0=TRUE, parallel=4)$simple,
      0)
    testthat::expect_identical(
      coverage(Ddflist_start1[[3]], t=5, isInteger=FALSE, start0=TRUE, parallel=4)$simple,
      0)
    testthat::expect_identical(
      coverage(Dlist_start1[[4]], t=6, isInteger=TRUE, start0=FALSE, parallel=1)$simple,
      0)
    testthat::expect_identical(
      coverage(Ddflist[[4]], t=5, isInteger=FALSE, start0=TRUE, parallel=1),
      perfect)
    testthat::expect_error(
      coverage(Dlist_noninteger[[4]], t=5, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
})
})
