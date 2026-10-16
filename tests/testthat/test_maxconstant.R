library(CAs)
library(testthat)

# Tests for input validation error messages in flexible.R functions

# ========================================
# markflex() function tests
# ========================================
Dpaley <- paleyCA(3,8,2) ## has a row of ones
Dcharacter <- matrix(as.character(Dpaley), nrow=12)
Dwithduplicates <- rbind(0,1,Dpaley)
Dmixed <- rbind(c(0,0,0,1,1,1,2,1,0),
                c(1,1,1,0,0,0,0,2,1),
                c(0,1,0,1,0,1,1,0,2))
Ddataframe <- as.data.frame(Dpaley)

test_that("maxconstant handling of inputs",
          {
    expect_error(maxconstant( Dmixed ),
                 "function maxconstant works for symmetric CAs only")
    expect_error(maxconstant( Ddataframe))
    expect_error(maxconstant( Dcharacter))
    expect_error(maxconstant(Dpaley, verbose=13))
    expect_error(maxconstant(Dpaley, one_is_enough = "TRUE"))
    expect_error(maxconstant(Dpaley, dupcheck = "TRUE"))
    expect_error(maxconstant(Dpaley, remove = "TRUE"))
          })
test_that("maxconstant correctly handles duplicate rows", {
  ## no duplicates
  expect_equal(dim(maxconstant(Dpaley)),c(12,8))
  expect_equal(dim(maxconstant(Dpaley, dupcheck=TRUE)),c(12,8))
  expect_equal(dim(maxconstant(Dwithduplicates,
                               dupcheck = TRUE)),c(13,8))
  expect_equal(dim(maxconstant(Dwithduplicates,
                               dupcheck = FALSE)),c(14,8))
})
test_that("maxconstant correctly handles NA values in clique rows",
          {
            expect_all_equal(unlist(coverage(
              maxconstant(markflex(paleyCA(3,8,2),2)[-c(7,11:12),] ),
              2)), 1)
          })
