library(CAs)
test_that("smolenBigFrame is handled correctly", {
  klasse <- class(CAs:::smolenBigFrame)
  expect_true(klasse=="NULL" || "data.frame" %in% klasse)
  if ("data.frame" %in% klasse){
    expect_true(eCAN(3, 12, 4, frozen=TRUE)$CAN >=
                  eCAN(3, 12, 4, frozen=FALSE)$CAN)
    expect_true(eCAK(5, 550, 2, frozen=TRUE)$CAK <=
                  eCAK(5, 550, 2, frozen=FALSE)$CAK)
  }

  if (klasse=="NULL"){
    ## should only happen without an internet connection
    expect_error(eCAN(3, 12, 4, frozen=FALSE))
    expect_error(eCAK(5, 550, 2, frozen=FALSE))
  }
})
