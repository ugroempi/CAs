library(CAs)

test_that("Ns_add1 provides direct and add1 sizes", {
  A <- Ns_add1(2, 8, 3, maxadd=2)
  expect_identical(rownames(A), c("direct", "add1", "add2"))
  ## one add1 step from k-1 columns adds v*(v-1) runs for t=2
  expect_equal(unlist(A["add1", "eCAN"]),
               unname(Ns(2, 7, 3, noderive=TRUE, nofuse=TRUE, noadd1=TRUE)["eCAN"]) + 6)
})

test_that("Ns has an add1 entry that can be switched off", {
  N <- Ns(2, 8, 3)
  expect_true(any(grepl(":add1$", names(N))))
  expect_false(any(grepl(":add1", names(Ns(2, 8, 3, noadd1=TRUE)))))
  expect_false(any(grepl(":add1", names(Ns(2, 8, 3, exclude="add1")))))
})

test_that("bestCA builds add1 designs of the announced size with ingredient", {
  N <- Ns(2, 8, 3, noderive=TRUE, nofuse=TRUE)
  entry <- N[grepl(":add1$", names(N))]
  code <- CAs:::labelToCode(names(entry), 2, 8, 3)
  D <- eval(parse(text=code))
  expect_equal(nrow(D), unname(entry))
  expect_equal(ncol(D), 8)
  expect_equal(prod(unlist(coverage(D, 2))), 1)
})

test_that("bestCA records the ingredient of an add1 design", {
  ex <- c("DWYER", "NIST", "CKRS", "scphfCA", "cphfCA", "CK_doublingCA",
          "RMScopiesCA", "fuse", "derive")
  D <- bestCA(3, 7, 3, exclude=ex)
  expect_identical(attr(D, "constr"), "CK_NRB:add1")
  expect_identical(attr(D, "ingredients")$source[c("constr", "t", "k", "v")],
                   list(constr="CK_NRB", t=3, k=6, v=3))
})
