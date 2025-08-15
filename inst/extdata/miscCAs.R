## miscellaneous arrays
## arrays contained as separate objects
## info collected in miscCAcat at the very top
##
miscCAcat <- rbind(
  data.frame(t=3, k=8, v=2, N=12, fns="ca12.2.8", nconst=2, PCAstatus=0, hasNA=FALSE, comment="created from paleyCA(3,23)"),
  data.frame(t=5, k=6, v=3, N=243, fns="DoE.base::L243.3.6-1", nconst=3, PCAstatus=0, hasNA=FALSE, comment=""),
  data.frame(t=5, k=12, v=3, N=729, fns="DoE.base::L729.3.12-1", nconst=3, PCAstatus=0, hasNA=FALSE, comment=""),
  data.frame(t=2, k=15, v=4, N=26, fns="ca26.4.15", nconst=1, PCAstatus=11, hasNA=FALSE, comment="SCA with k1=11 for up to 13 columns"),
  data.frame(t=3, k=6, v=4, N=64, fns="DoE.base::L64.4.6-1", nconst=4, PCAstatus=0, hasNA=FALSE, comment=""),
  data.frame(t=5, k=6, v=4, N=1024, fns="DoE.base::L1024.4.6-1", nconst=4, PCAstatus=0, hasNA=FALSE, comment=""),
  data.frame(t=2, k=3, v=6, N=36, fns="DoE.base::L36.2.8.6.3[,9:11]-1", nconst=6, PCAstatus=0, hasNA=FALSE, comment=""),
  data.frame(t=2, k=4, v=10, N=100, fns="DoE.base::L100.2.4.10.4[,5:8]-1", nconst=9, PCAstatus=3, hasNA=FALSE, comment=""),
  data.frame(t=2, k=7, v=12, N=144, fns="DoE.base::L144.12.7-1", nconst=9, PCAstatus=6, hasNA=FALSE, comment=""),
  data.frame(t=3, k=6, v=12, N=1728, fns="oa1728.12.6", nconst=1, PCAstatus=4, hasNA=FALSE, comment=""),
  data.frame(t=3, k=6, v=15, N=3375, fns="oa3375.15.6", nconst=1, PCAstatus=4, hasNA=FALSE, comment=""),
  data.frame(t=3, k=6, v=21, N=9261, fns="oa9261.21.6", nconst=1, PCAstatus=4, hasNA=FALSE, comment="")
)

## arrays from DoE.base can gain constant rows via maxconstant or
##      - where there are fewer than v constant rows - a PCA structure via CA_to_PCA
## PCAstatus=0 means that there are v constant rows - otherwise,
##      PCAstatus must be at least 1, as it is always possible to arrange rows such that
##      the first column has v distinct elements

## ca12.2.8
## make one row constant
D <- maxconstant(paleyCA(3,23))
rowSums(D)
## keep 8 columns with row 2 equal to 1
D <- D[,which(D[2,]==1)[5:12]]
ca12.2.8 <- postopNCK(D, 3, fixrows=2, innerRetry = 3, seed=17190)
attr(ca12.2.8, "rowOrder") <- NULL
attr(ca12.2.8, "seed") <- NULL
attr(ca12.2.8, "Call") <- NULL
attr(ca12.2.8, "origin") <- c("from paleyCA(3,23) by making one row constant,",
                              "picking last eight columns with second row also constant,",
                              "and reducing to 12 runs via postopNCK with seed 17190")

## ca26.4.15
hilf <- strsplit("220033101323300 030210112330122 113300232313302 003123221103122 132102303103211 101231020033333 202211323301202 203211032112311 312203001230220 321130332202130 212300123002123 332122313000332 330021223211110 303312110223010 231122010331131 321033210132201 021301003322012 01032013102⋆211 130122301222103 012032330211322 010300111121033 123013202011233 000000000000000 111111111110020 222222222220021 333333333330013", " ")
ca26.4.15 <- CAs:::funmakefromstrings(hilf[[1]])[c(23:26,1:22),]
     ## remark: PCA with k1=11 for up to 13 columns
sum(is.na(ca26.4.15)) ## one NA is as expected
table(ca26.4.15[,12])
## fix the CA (markflex easily recovers it)
ca26.4.15[is.na(ca26.4.15)] <- 1 ## avoid 0 for balance
coverage(ca26.4.15, 2)
class(ca26.4.15) <- c("ca", class(ca26.4.15))
attr(ca26.4.15, "t") <- 2
attr(ca26.4.15, "source") <- "simulated annealing (Torres-Jimenez)"
attr(ca26.4.15, "origin") <- "found in Colbourn and Torres-Jimenez (2013, Figure 1)"

################################################################################
## ca12.3.7 from postopNCK applied to Figure 2 of Colbourn and Torres-Jimenez (2013)
## commented out, because available via TJcat
##
# hilf <- CAs:::funmakefromstrings(c("2222001", "0022222", "1121201", "0120110",
#                    "2210202", "2102122", "1200020", "0211121",
#                    "1012100", "2001210", "0000001", "1111012",
#                    "2222222", "2222011"))
# dim(hilf)
# (ca12.3.7 <- postopNCK(hilf, 2))
# dim(ca12.3.7)     ## now optimal
# attr(ca12.3.7, "rowOrder") <- NULL
# eCAN(2,7,3)
# ## make into PCA
# ca12.3.7 <- maxconstant(ca12.3.7)
# ca12.3.7 <- CA_to_PCA(ca12.3.7, tryhard=TRUE)

## oa1728.12.6 #################################################################
## Lemma 3.5 of Ji and Yin 2012 (last row is the adder)
## (12.4.1)-DM over Z6xZ2
#00 00 00 00 00 00 00 00 00 00 00 00
#00 01 10 11 20 21 30 31 40 41 50 51
#00 10 01 40 30 51 20 50 41 21 31 11
#00 11 40 20 41 01 50 21 31 10 51 30
#00 21 10 11 40 41 30 51 20 01 50 31

## One has to use the separate elements in the construction, not Z12
## The paper has this construction for
## v=10,12,14,15,18,21,22, 24 (Table 2 with subsequent special cases 15 and 21)
## The Colbourn tables also have an OA for v=20, attributed to Ji and Yin, where???
## 12, 24, 15 and 21 are simple, the others require to construct an ingredient

## from Lemma 3.5 of Ji and Yin
DM12 <- array(NA, dim=c(4,12,2))
## the Z6 portion
DM12[1:4,1:12,1] <- rbind(
  rep(0,12),
  rep(0:5, each=2),
  c(0,1,0,4,3,5,2,5,4,2,3,1),
  c(0,1,4,2,4,0,5,2,3,1,5,3)
)
## the Z2 portion
DM12[1:4,1:12,2] <- rbind(
  rep(0,12),
  rep(0:1,6),
  c(0,0,1,0,0,1,0,0,1,1,1,1),
  c(0,1,0,0,1,1,0,1,1,0,1,0)
)
dimnames(DM12) <- list(NULL, NULL, group=c("Z6","Z2"))
#00 21 10 11 40 41 30 51 20 01 50 31
adder12 <- list(Z6=c(0,2,1,1,4,4,3,5,2,0,5,3), Z2=rep(0:1,6))

## from the proof of Theorem 3.2 of Ji and Yin
# (d1j + u, d2j + u, d3j + u + e + sj , d4j + u + e + sj , e, e + sj )

## this implementation is very adhoc for the single file
## should be made into a general functionality, then I can also produce the OA for six 24-level
##    columns which has 13824 runs, or any other construction based on Zm x Zn fields
additionTable <-
  rbind(cbind(circular(0:5), circular(6:11)), cbind(circular(6:11), circular(0:5)))
additionTable <- additionTable[,c(1,6:2,7,12:8)]

multiplicationTable <- cbind(rep(0,12),
                             c(0:5,0:5),
                             rep(c(0,2,4),4),
                             rep(c(0,3),6),
                             rep(c(0,4,2), 4),
                             rep(c(0,5:1), 2),
                             rep(c(0,6), each=6),
                             0:11,
                             c(0,2,4,0,2,4,6,8,10,6,8,10),
                             c(0,3,0,3,0,3,6,9,6,9,6,9),
                             c(0,4,2,0,4,2,6,10,8,6,10,8),
                             c(0,5:1, 6, 11:7))

plus12 <- function(x, y) additionTable[x%%12 + 1, y%%12 + 1]
mal12 <- function(x,y) multiplicationTable[x%%12 + 1, y%%12 + 1]

addernum12 <- 6*adder12$Z2 + adder12$Z6

## corresponds to the mapping 6*Z2 + Z6

n <- 12
nc <- 6
Gs <- c(6,2)

oa <- matrix(NA, nrow = n^3, ncol = nc)
row_idx <- 1

for (u in 0:(n-1)) {
  for (e in 0:(n-1)) {
    for (j in 1:ncol(DM12)){
      nowcol <- c(DM12[,j,]%*%c(1,6)) ## make numeric
      nows <- addernum12[j] # sapply(adder12, function(obj) obj[j])
      # Use col j of the DCA and element j of adder
        oa[row_idx, ] <- c(
          plus12(nowcol[1], u),
          plus12(nowcol[2], u),
          plus12(plus12(plus12(nowcol[3], u), e), nows),
          plus12(plus12(plus12(nowcol[4], u), e), nows),
          e,
          plus12(e, nows)
        )
      row_idx <- row_idx + 1
    }
  }
}

dim(oa)  # Should be 1728 x 6

# DoE.base::GWLP(oa, kmax=4)  # 1   0   0   0 165

oa1728.12.6 <- maxconstant(oa, one_is_enough = TRUE)
oa1728.12.6 <- oa1728.12.6[,c(2:4,6,1,5)]
## manually bring to PCA shape, as CA_to_PCA takes very long in spite of the obvious case
for (i in 1:4) oa1728.12.6 <- permvals(oa1728.12.6, i, oa1728.12.6[1:12,i], 0:11)
class(oa1728.12.6) <- c("ca", "oa", class(oa1728.12.6))
attr(oa1728.12.6, "t") <- 3
attr(oa1728.12.6, "origin") <- "Lemma 3.5 of Ji and Yin 2012"
attr(oa1728.12.6, "source") <- "orthogonal array (Ji-Yin)"
attr(oa1728.12.6, "PCAstatus") <- is.PCA(oa1728.12.6)

##########################################################################

## oa .15.6 ##############################################################
## Lemma 4.1
## (15, 4, 1) DM over Z_15
DM15 <- rbind(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
c(0, 2, 7, 1, 11, 4, 10, 13, 3, 6, 12, 14, 5, 9, 8),
c(0, 10, 9, 8, 1, 7, 4, 2, 14, 5, 13, 12, 11, 6, 3))

addernum15 <- c(0, 9, 1, 4, 2, 14, 7, 12, 6, 8, 10, 5, 11, 3, 13)

n <- 15
nc <- 6
Gs <- n

oa <- matrix(NA, nrow = n^3, ncol = nc)
row_idx <- 1

for (u in 0:(n-1)) {
  for (e in 0:(n-1)) {
    for (j in 1:ncol(DM15)){
      nowcol <- DM15[,j] ## make numeric
      nows <- addernum15[j] # sapply(adder15, function(obj) obj[j])
      # Use col j of the DCA and element j of adder
      oa[row_idx, ] <- c(
        (nowcol[1] + u)%%15,
        (nowcol[2]+ u)%%15,
        (((nowcol[3]+ u)+ e)+ nows)%%15,
        (((nowcol[4] + u) + e) + nows)%%15,
        e,
        (e + nows)%%15
      )
      row_idx <- row_idx + 1
    }
  }
}

dim(oa)  # Should be 3375 x 6
coverage(oa,3)
oa <- oa[,c(2:4,6,1,5)]
for (i in 1:4) oa <- CAs:::permvals(oa, i, oa[1:15,i], 0:14)
head(oa, n=16)
oa3375.15.6 <- oa
class(oa3375.15.6) <- c("ca", class(oa3375.15.6))
attr(oa3375.15.6, "source") <- "orthogonal array Ji-Yin"
attr(oa3375.15.6, "origin") <- "Ji and Yin Lemma 4.1"
attr(oa3375.15.6, "t") <- 3
attr(oa3375.15.6, "PCAstatus") <- list(type="SCA", k1=4, k2=2)

DM21 <- rbind(
c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    c(2, 4, 6, 12, 14, 11, 18, 1, 13, 16),
      c(3, 16, 13, 17, 1, 7, 12, 11, 15, 19)
)
DM21 <- cbind(DM21, (-DM21)%%21, 0)
dim(DM21)
### the adder in the paper has only 20 elements
### how to fix?
### presumably, a 0 must be inserted somewhere
### neither the end nor the start works
addernum21 <- c(5, 11, 13, 14, 3, 2, 6, 9, 20, 4)
addernum21 <- c(addernum21, (-addernum21)%%21, 0)

n <- 21
nc <- 6

oa <- matrix(NA, nrow = n^3, ncol = nc)
row_idx <- 1

for (u in 0:(n-1)) {
  for (e in 0:(n-1)) {
    for (j in 1:ncol(DM21)){
      nowcol <- DM21[,j]
      nows <- addernum21[j]
      # Use col j of the DCA and element j of adder
      oa[row_idx, ] <- c(
        (nowcol[1] + u)%%21,
        (nowcol[2] + u)%%21,
        (((nowcol[3]+ u)+ e)+ nows)%%21,
        (((nowcol[4] + u) + e) + nows)%%21,
        e,
        (e + nows)%%21
      )
      row_idx <- row_idx + 1
    }
  }
}

dim(oa)  # Should be 9261 x 6
coverage(oa,3)
oa <- oa[,c(2:4,6,1,5)]
for (i in 1:4) oa <- CAs:::permvals(oa, i, oa[1:21,i], 0:20)
head(oa, n=22)
oa9261.21.6 <- oa
class(oa9261.21.6) <- c("ca", class(oa9261.21.6))
attr(oa9261.21.6, "source") <- "orthogonal array Ji-Yin"
attr(oa9261.21.6, "origin") <- "Ji and Yin Lemma 4.1"
attr(oa9261.21.6, "t") <- 3
attr(oa9261.21.6, "PCAstatus") <- list(type="SCA", k1=4, k2=2)

# save(ca12.2.8, file="d:/rtests/CAs/data/ca12.2.8.rda", compress="xz")
# save(ca26.4.15, file="d:/rtests/CAs/data/ca26.4.15.rda", compress="xz")
# save(oa1728.12.6, file="d:/rtests/CAs/data/oa1728.12.6.rda", compress="xz")
# save(oa3375.15.6, file="d:/rtests/CAs/data/oa3375.15.6.rda", compress="xz")
# save(oa9261.21.6, file="d:/rtests/CAs/data/oa9261.21.6.rda", compress="xz")
# save(miscCAcat, file="d:/rtests/CAs/data/miscCAcat.rda", compress="xz")
