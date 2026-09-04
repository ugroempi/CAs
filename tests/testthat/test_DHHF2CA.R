library(CAs)
test_that("DHHF2CA works", {
  ## check for wrong number of arrays, did not manage to check for the error message
  expect_error(DHHF2CA(t(SCA_Busht(3,3))[1:3,],list(bestCA(3,27,2), bestCA(3,24,2))))
  ## create an array with chi=0
    myDHHF <- t(SCA_Busht(11,2))[1:4,]
    myDHHF <- myDHHF[,which(myDHHF[4,]<8)] ## last row 8 only
    Turan(3,4) ## needs more rows, i.e., 4
    D1 <- compositCA(3,11,4) ## has a single constant row (level 3)
    D2 <- CK_NRB(8,4)        ## has four constant rows
    aus <- DHHF2CA(myDHHF, list(D1,D2))
    hilf <- caverify::ca_verify(aus,3)
    expect_true(hilf$covered)
    expect_equal(hilf$N, 513)
    ## create an array with chi=1
    ## reduce two further rows to 8 levels
    myDHHF <- myDHHF[,which(myDHHF[3,]<8)]
    myDHHF <- myDHHF[,which(myDHHF[2,]<8)]
    aus <- DHHF2CA(myDHHF, list(D1,D2))
    hilf <- caverify::ca_verify(aus,3)
    expect_true(hilf$covered)
    expect_equal(hilf$N, 396)
    ## three blocks tested with powerCA
})
