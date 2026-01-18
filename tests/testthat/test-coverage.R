library(CAs)
test_that("coverage functionality works", {
         perfect <- coverage(bestCA(3,8,4),3)
         covNA <- perfect
         for (i in 1:4) covNA[[i]] <- NA

    ## a simple mixed-level CA of strength 2
    ## that has also a few flexible values
    D <- MCA2(nlevels=c(6,4,3,2,2,2,2,2,2), seed=2481) ## 24 runs
    cD1 <- coverage(D, 2, verbose=1)
    cD2 <- coverage(D, 2, verbose=2)
    cD13 <- coverage(D, 3, verbose=1)
    cD23 <- coverage(D, 3, verbose=2)

    expect_s3_class(cD1, "CAcoverage")
    expect_s3_class(cD2, "CAcoverage")
    expect_s3_class(cD13, "CAcoverage")
    expect_s3_class(cD23, "CAcoverage")

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
    expect_type(cD23$tabs, "list")

    ## check coverplot output
    expect_error(coverplot(cD13),
                 "If t is not given, coverageD must be a coverage object from function coverage with verbose=2")
    expect_identical(length(coverplot(cD13, t=3, plot=FALSE)), 4L)
    expect_error(coverplot(coverage(D,3)),
          "If t is not given, coverageD must be a coverage object from function coverage with verbose=2")
    expect_error(coverplot(cD13),
          "If t is not given, coverageD must be a coverage object from function coverage with verbose=2")
    expect_error(coverplot(coverage(D,3), t=3),
          "coverageD must be a coverage object from function coverage with verbose at least 1")
    expect_error(coverplot(D), "coverageD must have class 'CAcoverage'. Did you want to specify a design D?")
    expect_error(coverplot(D=D), "If coverageD is missing, both D and t must be specified.")
    expect_true(length(coverplot(D=D, t=3, plot=FALSE))==4L)
    expect_identical(names(coverplot(D=D, t=3, plot=FALSE)),
                     c("x","y","t","coverage"))
    expect_identical(names(coverplot(cD23, plot=FALSE)),
                     c("x","y","t","coverage"))
    expect_equal(coverplot(D=D, t=3, plot=FALSE),
                     coverplot(cD23, plot=FALSE))
    ## with default plot=TRUE
    # expect_equal(coverplot(D=D, t=3), coverplot(cD23))
    # expect_equal(coverplot(D=D, t=3, type="tuples"),
    #             coverplot(cD23, type="tuples"))

    #### now for coverage_iter
    cD1perfect <- coverage_iter(D, 2, abortnot1=TRUE)
    cD2perfect <- coverage_iter(D, 2, abortnot1=FALSE)
    cD1not1 <- coverage_iter(D, 3, abortnot1=TRUE)
    cD2not1 <- coverage_iter(D, 3, abortnot1=FALSE)

    expect_s3_class(cD1perfect, "CAcoverage")
    expect_s3_class(cD2perfect, "CAcoverage")
    expect_s3_class(cD1not1, "CAcoverage")
    expect_s3_class(cD2not1, "CAcoverage")

    expect_identical(length(cD1perfect), 4L)
    expect_identical(length(cD2perfect), 4L)
    expect_identical(length(cD1not1), 5L)
    expect_identical(length(cD2not1), 4L)

    ##
    expect_identical(cD1perfect[1:4], unclass(perfect))
    expect_identical(cD2perfect[1:4], unclass(perfect))
    expect_identical(cD1not1[1:4], unclass(covNA))
    expect_equal(cD2not1$total*332, 179)
    expect_identical(cD2not1$ave*1728, 1103)
    expect_identical(cD2not1$min*24, 7)
    expect_identical(cD2not1$simple*21, 1)

    ### with data.frame D
    D <- as.data.frame(D)
    cD1 <- coverage(D, 2, verbose=1)
    cD2 <- coverage(D, 2, verbose=2)
    cD13 <- coverage(D, 3, verbose=1)
    cD23 <- coverage(D, 3, verbose=2)

    expect_s3_class(cD1, "CAcoverage")
    expect_s3_class(cD2, "CAcoverage")
    expect_s3_class(cD13, "CAcoverage")
    expect_s3_class(cD23, "CAcoverage")

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
    expect_type(cD23$tabs, "list")

    ## now for coverage_iter
    cD1perfect <- coverage_iter(D, 2, abortnot1=TRUE)
    cD2perfect <- coverage_iter(D, 2, abortnot1=FALSE)
    cD1not1 <- coverage_iter(D, 3, abortnot1=TRUE)
    cD2not1 <- coverage_iter(D, 3, abortnot1=FALSE)

    expect_s3_class(cD1perfect, "CAcoverage")
    expect_s3_class(cD2perfect, "CAcoverage")
    expect_s3_class(cD1not1, "CAcoverage")
    expect_s3_class(cD2not1, "CAcoverage")

    expect_identical(length(cD1perfect), 4L)
    expect_identical(length(cD2perfect), 4L)
    expect_identical(length(cD1not1), 5L)
    expect_identical(length(cD2not1), 4L)

    ##
    expect_identical(cD1perfect[1:4], unclass(perfect))
    expect_identical(cD2perfect[1:4], unclass(perfect))
    expect_identical(cD1not1[1:4], unclass(covNA))
    expect_equal(cD2not1$total*332, 179)
    expect_identical(cD2not1$ave*1728, 1103)
    expect_identical(cD2not1$min*24, 7)
    expect_identical(cD2not1$simple*21, 1)

    ## tabs as a matrix for uniform CAs with verbose=2
    cD <- coverage(KSK(k=12), 3, verbose=2)
    expect_true(length(cD$ncovereds)==220L)
    expect_type(cD$tabs, "integer")
    expect_true(is.matrix(cD$tabs))
    expect_true(ncol(cD$tabs)==220L)
    expect_true(nrow(cD$tabs)==8L) ## for the combinations

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
##    factor H: parallel=1 or parallel=2
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

    test_that("CA-based test for coverage inputs (uniform inputs only)", {
    expect_error(
      coverage(Ddflist_start1[[4]], t=6, isInteger=TRUE, start0=TRUE, parallel=2),
      regexp="start0 is TRUE, but columns do not start at 0"
    )
    expect_identical(
      coverage(Dlist[[2]], t=4, isInteger=FALSE, start0=FALSE, parallel=2)$simple,
      0)
    expect_identical(
      coverage(Dlist_start1[[1]], t=2, isInteger=TRUE, start0=FALSE, parallel=2),
      perfect)
    expect_identical(
      coverage(Ddflist_noninteger[[4]], t=6, isInteger=FALSE, start0=FALSE, parallel=2)$simple,
      0)
    expect_error(
      coverage(Dlist_noninteger[[3]], t=5, isInteger=TRUE, start0=TRUE, parallel=2),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    expect_identical(
      coverage(Dlist_noninteger[[2]], t=4, isInteger=FALSE, start0=TRUE, parallel=1)$simple,
      0)
    expect_identical(
      coverage(Dlist_noninteger[[4]], t=5, isInteger=FALSE, start0=TRUE, parallel=2),
      perfect)
    expect_identical(
      coverage(Ddflist[[2]], t=3, isInteger=TRUE, start0=TRUE, parallel=2),
      perfect)
## modified isInteger setting for the following row,
## in order to have also error with start0=FALSE
    expect_error(
      coverage(Ddflist[[3]], t=4, isInteger=TRUE, start0=FALSE, parallel=2),
      regexp="start0 is FALSE, but columns do not start at 1"
    )
    expect_error(
      coverage(Dlist_start1[[3]], t=4, isInteger=TRUE, start0=TRUE, parallel=1),
      regexp="start0 is TRUE, but columns do not start at 0"
    )
    expect_error(
      coverage(Ddflist_noninteger[[1]], t=2, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    expect_identical(
      coverage(Dlist[[1]], t=3, isInteger=FALSE, start0=FALSE, parallel=1)$simple,
      0)
    expect_identical(
      coverage(Ddflist[[1]], t=3, isInteger=TRUE, start0=TRUE, parallel=1)$simple,
      0)
    expect_identical(
      coverage(Dlist_noninteger[[3]], t=4, isInteger=FALSE, start0=FALSE, parallel=1),
      perfect)
    expect_identical(
      coverage(Ddflist_noninteger[[1]], t=2, isInteger=FALSE, start0=TRUE, parallel=2),
      perfect)
    expect_error(
      coverage(Ddflist_noninteger[[2]], t=4, isInteger=TRUE, start0=TRUE, parallel=2),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    expect_identical(
      coverage(Ddflist_start1[[2]], t=3, isInteger=FALSE, start0=FALSE, parallel=1),
      perfect)
    expect_error(
      coverage(Dlist_noninteger[[2]], t=3, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    expect_error(
      coverage(Ddflist_noninteger[[3]], t=5, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
    expect_identical(
      coverage(Dlist_noninteger[[1]], t=3, isInteger=FALSE, start0=TRUE, parallel=2)$simple,
      0)
    expect_identical(
      coverage(Ddflist_start1[[3]], t=5, isInteger=FALSE, start0=TRUE, parallel=2)$simple,
      0)
    expect_identical(
      coverage(Dlist_start1[[4]], t=6, isInteger=TRUE, start0=FALSE, parallel=1)$simple,
      0)
    expect_identical(
      coverage(Ddflist[[4]], t=5, isInteger=FALSE, start0=TRUE, parallel=1),
      perfect)
    expect_error(
      coverage(Dlist_noninteger[[4]], t=5, isInteger=TRUE, start0=FALSE, parallel=1),
      regexp="isInteger is TRUE, but D has non-numeric content"
    )
})

### now for coverage_iter with factor H: abortnot1 instead of parallel
  expect_error(
    coverage_iter(Ddflist_start1[[4]], t=6, isInteger=TRUE, start0=TRUE, abortnot1=TRUE),
    regexp="start0 is TRUE, but columns do not start at 0"
  )
  suppressMessages(expect_identical(
    coverage_iter(Dlist[[2]], t=4, isInteger=FALSE, start0=FALSE, abortnot1=TRUE)$message,
    "First projection with strength 4 coverage violated: 1:2:3:4"))
  expect_identical(
    coverage_iter(Dlist_start1[[1]], t=2, isInteger=TRUE, start0=FALSE, abortnot1=TRUE),
    perfect)
  suppressMessages({
    hilf <- coverage_iter(Ddflist_noninteger[[4]], t=6, isInteger=FALSE, start0=FALSE, abortnot1=TRUE)
    expect_identical(hilf$message,
                     "First projection with strength 6 coverage violated: 1:2:3:4:5:6")
    expect_identical(unclass(hilf[1:4]), unclass(covNA))
  })
  expect_error(
    coverage_iter(Dlist_noninteger[[3]], t=5, isInteger=TRUE, start0=TRUE, abortnot1=TRUE),
    regexp="isInteger is TRUE, but D has non-numeric content"
  )
  expect_identical(
    coverage_iter(Dlist_noninteger[[2]], t=4, isInteger=FALSE, start0=TRUE, abortnot1=FALSE)$simple,
    0)
  expect_identical(
    coverage_iter(Dlist_noninteger[[4]], t=5, isInteger=FALSE, start0=TRUE, abortnot1=TRUE),
    perfect)
  expect_identical(
    coverage_iter(Ddflist[[2]], t=3, isInteger=TRUE, start0=TRUE, abortnot1=TRUE),
    perfect)
  ## modified isInteger setting for the following row,
  ## in order to have also error with start0=FALSE
  expect_error(
    coverage_iter(Ddflist[[3]], t=4, isInteger=TRUE, start0=FALSE, abortnot1=TRUE),
    regexp="start0 is FALSE, but columns do not start at 1"
  )
  expect_error(
    coverage_iter(Dlist_start1[[3]], t=4, isInteger=TRUE, start0=TRUE, abortnot1=FALSE),
    regexp="start0 is TRUE, but columns do not start at 0"
  )
  expect_error(
    coverage_iter(Ddflist_noninteger[[1]], t=2, isInteger=TRUE, start0=FALSE, abortnot1=FALSE),
    regexp="isInteger is TRUE, but D has non-numeric content"
  )
  expect_identical(
    coverage_iter(Dlist[[1]], t=3, isInteger=FALSE, start0=FALSE, abortnot1=FALSE)$simple,
    0)
  expect_identical(
    coverage_iter(Ddflist[[1]], t=3, isInteger=TRUE, start0=TRUE, abortnot1=FALSE)$simple,
    0)
  expect_identical(
    coverage_iter(Dlist_noninteger[[3]], t=4, isInteger=FALSE, start0=FALSE, abortnot1=FALSE),
    perfect)
  expect_identical(
    coverage_iter(Ddflist_noninteger[[1]], t=2, isInteger=FALSE, start0=TRUE, abortnot1=TRUE),
    perfect)
  expect_error(
    coverage_iter(Ddflist_noninteger[[2]], t=4, isInteger=TRUE, start0=TRUE, abortnot1=TRUE),
    regexp="isInteger is TRUE, but D has non-numeric content"
  )
  expect_identical(
    coverage_iter(Ddflist_start1[[2]], t=3, isInteger=FALSE, start0=FALSE, abortnot1=FALSE),
    perfect)
  expect_error(
    coverage_iter(Dlist_noninteger[[2]], t=3, isInteger=TRUE, start0=FALSE, abortnot1=FALSE),
    regexp="isInteger is TRUE, but D has non-numeric content"
  )
  expect_error(
    coverage_iter(Ddflist_noninteger[[3]], t=5, isInteger=TRUE, start0=FALSE, abortnot1=FALSE),
    regexp="isInteger is TRUE, but D has non-numeric content"
  )
  suppressMessages(expect_identical(
    coverage_iter(Dlist_noninteger[[1]], t=3, isInteger=FALSE, start0=TRUE, abortnot1=TRUE)$message,
    "First projection with strength 3 coverage violated: 1:2:3"))
  suppressMessages(expect_identical(
    coverage_iter(Ddflist_start1[[3]], t=5, isInteger=FALSE, start0=TRUE, abortnot1=TRUE)$message,
    "First projection with strength 5 coverage violated: 1:2:3:4:5"))
  expect_identical(
    coverage_iter(Dlist_start1[[4]], t=6, isInteger=TRUE, start0=FALSE, abortnot1=FALSE)$simple,
    0)
  expect_identical(
    coverage_iter(Ddflist[[4]], t=5, isInteger=FALSE, start0=TRUE, abortnot1=FALSE),
    perfect)
  expect_error(
    coverage_iter(Dlist_noninteger[[4]], t=5, isInteger=TRUE, start0=FALSE, abortnot1=FALSE),
    regexp="isInteger is TRUE, but D has non-numeric content"
  )
})
