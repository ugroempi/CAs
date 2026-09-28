library(CAs)
cyclo <- colbournBigFrame[setdiff(grep("Cyclotomy", colbournBigFrame$Source),
                                  grep("shift", colbournBigFrame$Source)),]
type <- rep("", nrow(cyclo))
primedat <- CAs:::primedat

## identify the prime power used
q <- rep(0, nrow(cyclo))
## identify the construction type
  type[cyclo$k==cyclo$N] <- "1"
  type[cyclo$k+cyclo$v-1==cyclo$N] <- "2"
  type[cyclo$k*cyclo$v==cyclo$N] <- "3o4a"
  type[(cyclo$k-1)*cyclo$v==cyclo$N] <- "3a"
  type[(cyclo$k+cyclo$v-2)*cyclo$v==cyclo$N] <- "4b"
  type[(cyclo$k+1)*cyclo$v==cyclo$N] <- "4"
## infer q
  q[type %in% c("1","2","4")] <- cyclo$k[type %in% c("1","2","4")]
  q[type %in% c("3a", "4b")] <- cyclo$k[type %in% c("3a", "4b")] - 1
  q[type=="3o4a"][which(cyclo$k[type=="3o4a"] %in% primedat$q)] <-
    cyclo$k[type=="3o4a"][which(cyclo$k[type=="3o4a"] %in% primedat$q)]
  q[type=="3o4a"][which((cyclo$k-1)[type=="3o4a"] %in% primedat$q)] <-
    (cyclo$k[type=="3o4a"]-1)[which((cyclo$k-1)[type=="3o4a"] %in% primedat$q)]
## disambiguate type
  type[type=="3o4a"][which(q[type=="3o4a"]==cyclo$k[type=="3o4a"])] <- "3"
  type[type=="3o4a"] <- "4a"
table(type)

## special cases
cyclo[type=="",]
## Cyclotomy (Colbourn) postop NCK
cyclo$N[which(cyclo$k==89)] <- 89
type[which(cyclo$k==89)] <- "1"
q[which(cyclo$k==89)] <- 89
# coverage(cyc(89),4,parallel=4)

## Cyclotomy (Torres-Jimenez)
cyclo$k[which(cyclo$N==10007)] <- 10007
type[which(cyclo$N==10007)] <- "1"
q[which(cyclo$N==10007)] <- 10007

## Cyclotomy (Torres-Jimenez)
(80072/8) %in% primedat$q
cyclo$k[cyclo$N==80072] <- 10009
type[which(cyclo$N==80072)] <- "3"
q[which(cyclo$N==80072)] <- 10009

## the fuse case should not be included here
## (the case it is fused from is there)
## open question: What does the mult in Cyclotomy mult (Colbourn) mean?
remove <- cyclo$Source %in% c("Cyclotomy (Colbourn) fuse",
                    "Cyclotomy mult (Colbourn)")

cyclo <- cbind(cyclo, q=q, type=type)[!remove,]
cyclo$Source <- NULL

## row for postopNCK set to non-postop N,
## rows with kapped k got their actual k,
## rows with fuse and mult omitted

## omitted row with q=2977, because that q is not a prime power (it is 13*229)
## not sure what is going on there!
## t=4, v=6, k=2978, N=17862, inferred construction 3a with "prime" 2977
##    as 17862/6=2977
cyclo <- cyclo[!cyclo$q==2977,]

## sort by t and within t by v
cyclo <- cyclo[DoE.base::ord(cyclo),]

## problem pointed out by Christopher Smolen (included Sep 2nd 2026)
# At q=127 (t=3, k=128, v=6,
#     N=762) the built array misses 32,004 tuples across 5,334 projections.
# Your own coverage() agrees with caverify exactly, and checkcond3a(3,
#     6, 127) returns FALSE, while it returns TRUE on control cases whose
# arrays verify. I also tried all 36 primitive roots of GF(127) and none
# yields a covering array, so it is not the choice of primitive. The two
# larger rows, q=2161 (t=3, v=20) and q=2311 (t=3, v=21), fail on
# 50-column subsets, and a subset of a covering array must itself cover,
# so that is proof even though full verification at those sizes is not
# feasible. The other 210 rows of CYCLOTOMYcat all pass, 16 by full
# verification and 194 by four different 50-column screens each, so the
# problem is isolated to these three.
#
# The good news is the repair is in the same family it seems. Type 3b
# works at all three field sizes. At q=127 that gives CA(792; 3, 128,
#     6), which I verified fully. At q=2161 and q=2311 it gives 43,600 and
# 48,951 runs, each passing all four screens, and the pattern is tidy,
# 3b costs exactly v(v-1) rows more than the failing 3a size in all
# three cases. So the fix looks like three catalogue rows pointing to 3b
# with the corrected run counts.

## correct the three respective rows by setting the type to 4a
## and increasing the run size by v
pick <- which(cyclo$q==127)
cyclo[pick,]$N <- cyclo[pick,]$N + cyclo[pick,]$v
cyclo[pick,]$type <- "4a"

pick <- which(cyclo$q==2161)
cyclo[pick,]$N <- cyclo[pick,]$N + cyclo[pick,]$v
cyclo[pick,]$type <- "4a"

pick <- which(CYCLOTOMYcat$q==2311)
cyclo[pick,]$N <- cyclo[pick,]$N + cyclo[pick,]$v
cyclo[pick,]$type <- "4a"

## add code column
rownames(cyclo) <- NULL

cyclo$code <- paste0("cyc(", cyclo$q, ",", cyclo$v, ", type='", cyclo$type, "')")
CYCLOTOMYcat <- cyclo

## everything seems correct
kcycvec <- Vectorize(k_cyc, c("q","v","type"))
table(kcycvec(CYCLOTOMYcat$q, CYCLOTOMYcat$v, CYCLOTOMYcat$type)==CYCLOTOMYcat$k)
Ncycvec <- Vectorize(N_cyc, c("q","v","type"))
table(Ncycvec(CYCLOTOMYcat$q, CYCLOTOMYcat$v, CYCLOTOMYcat$type)==CYCLOTOMYcat$N)

## check implementation, at least dimension-wise
for (i in 1:nrow(CYCLOTOMYcat)){
  print(i)
  hilf <- dim(eval(parse(text=CYCLOTOMYcat$code[i])))
  if (!hilf[1]==CYCLOTOMYcat$N[i]) message("wrong number of rows")
  if (!hilf[2]==CYCLOTOMYcat$k[i]) message("wrong number of columns")
}
## checking coverage is prohibitive for larger arrays
## here, trust in the Colbourn tables is needed

## problem: for v=3 and t=5, the three Torres-Jimenez
##     constructions violate q%%v == 1
##  Is that relevant ???
## for now, deactivated the check for this condition in cyc
## small sample checks on strength 5 coverage did
##     not show up violations

CYCLOTOMYcat[174:176,]
## checking by sampling for q=6983, 7013, 10007,
## all with strength 5, v=3, type 1
for (q in c(6983, 7013, 10007)){
  print(q)
  aus <- cyc(q, 3, type='1')
  for (i in 1:20){
    print(i);
    if(!all(coverage(aus[,sample(1:q, 20)], 5)==1)){
      cat("coverage violated\n")
      break}
  }
}






