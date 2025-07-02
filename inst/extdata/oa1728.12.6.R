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

DoE.base::GWLP(oa, kmax=4)
oa1728.12.6 <- oa

