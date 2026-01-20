library(CAs)
test_that("cyclotomy functionality works", {
  ## use all types
  ## use 2 levels with type 4b, more than 2 levels with the other types
  ## figure out the rows of CYCLOTOMYcat to use
  which.types <- lapply(names(table(CYCLOTOMYcat$type)),
                    function(obj) which(CYCLOTOMYcat$type==obj))
  which.typesvg2 <- lapply(which.types[1:6], function(obj) obj[which(CYCLOTOMYcat$v[obj]>2)])
  which.types[1:6] <- which.typesvg2
  which.smallest <- sapply(which.types,
                           function(obj) obj[which.min(CYCLOTOMYcat$N[obj])])
  names(which.smallest) <- names(table(CYCLOTOMYcat$type))
  (N.smallest <- CYCLOTOMYcat$N[which.smallest])
  (t.smallest <- CYCLOTOMYcat$t[which.smallest])
  (k.smallest <- CYCLOTOMYcat$k[which.smallest])
  (v.smallest <- CYCLOTOMYcat$v[which.smallest])
  ## create a design with each of 7 types
  ## (the smallest possible for v>2 in all cases except 4b)
  ## with the relevant t of that smallest design
  ## and the maximum k or k one smaller or k larger than maximum possible for the t-v combination
k.toolarge <- sapply(1:7, function(obj)
  max(CYCLOTOMYcat$k[which(CYCLOTOMYcat$t>=t.smallest[obj] &
                       CYCLOTOMYcat$v==v.smallest[obj])])) + 1
deslist <- vector(mode="list")
deserrlist <- vector(mode="list")
for (i in 1:7){
  deslist[[i]] <- cyclotomyCA(t.smallest[i], k.smallest[i], v.smallest[i])
  deserrlist[[i]] <- try(cyclotomyCA(t.smallest[i], k.toolarge[i],
                                     v.smallest[i]), silent = TRUE)
}
expect_equal(sapply(deslist, nrow), N.smallest)
expect_equal(sapply(deslist, ncol), k.smallest)
expect_equal(sapply(deslist, max), v.smallest - 1)
## all errors are the same
expect_equal(length(unique(deserrlist)), 1L)
## error message is as expected
expect_equal(deserrlist[[1]][1], "Error in cyclotomyCA(t.smallest[i], k.toolarge[i], v.smallest[i]) : \n  k is too large for the combination of t and v\n")
## simple errors
expect_error(cyclotomyCA(7,12,4), "t is too large")
expect_error(cyclotomyCA(5,12,30), "v is too large")
## non-error
expect_equal(cyclotomyCA(-3, 12,4), cyclotomyCA(2, 12,4))

#### also do v=8 (because of gf)
  # 3 8   194  1552   193   4a  cyc(193,8, type='4a')
  D <- cyclotomyCA(3, 194, 8)
  expect_equal(dim(D), c(1552, 194))
  expect_equal(max(D), 7)
})
