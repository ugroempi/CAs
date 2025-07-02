## cover starters Colbourn and Keri
## 2-level, strengths 4, 5, 6

######################### strength 4 ###############################
## lemma 2, k+1 columns in 2k rows
# 23 00000111101011001100101
# 27 000110100110111001110101000
# 28 0001101000100011110110011010
# 29 00011010001000101111001101011
# 30 000110100010001110111100101101
# 31 0001101000100011100111111010010
# 32 00011010001000110111110100101110
# 33 000110100010001101101011001111101
## lemma 1, k columns in 2k rows
# 21 000011101000100100001
# 24 000110011110010010101000
# 25 0001101010110010000000111
# 26 00011010000110101011100100

## lemma 2, k+1 columns in 2k rows
CS_CKStarters2 <- lapply(
  ## label is k, length is k-1
  list(
    '24'="00000111101011001100101",   ## 46 runs
    '28'="000110100110111001110101000",  ## 54 runs
    '29'= "0001101000100011110110011010",  ## 56 runs
    '30'= "00011010001000101111001101011",   ## 58 runs
    '31'= "000110100010001110111100101101",    ## 60 runs
    '32'= "0001101000100011100111111010010",
    '33'= "00011010001000110111110100101110",
    '34'= "000110100010001101101011001111101"),
  function(obj) as.numeric(unlist(strsplit(obj, "", fixed=TRUE)))
)

CS_CKStarters1 <- lapply(
  list(
    '21'= "000011101000100100001",  ## 42 runw
    '24'= "000110011110010010101000",  ## 48 runs
    '25'= "0001101010110010000000111",  ## 50 runs
    '26'= "00011010000110101011100100"),  ## 52 runs
  function(obj) as.numeric(unlist(strsplit(obj, "", fixed=TRUE)))
)

CS_CKStarters3 <- lapply(
  ## for Corollary 3, two additional constant rows needed
  ## versus lemma 2 construction
  ## only k=12 brings improvement, other cases covered with fewer runs
  ## for the other starters
  list(
    '12'= "00010110111",        ## 24 runs
    '20'= "0000101011110010011",  ## 40 runs
    '22'= "000100011111100101101",  ## 44 runs
    '23'= "0001110111011010011100",   ## 46 runs  compare to 23 columns of 24 of 2
                                              ##  2 has better 5-coverage
    '25'= "000111011101100010010110",   ## 50 runs  compare to 1, same
    '26'= "0001110011011011010101000",    ## 52 runs  compare to 1, same
    '27'= "00011101101001101110001010"),    ## 54 runs  compare to 27 columns of 28 of 2
                                               ##  3 has better 5-coverage
    function(obj) as.numeric(unlist(strsplit(obj, "", fixed=TRUE)))
  )

## arrange the best starters
## giving priority to 3 over the others for tied strength t+1 coverage
##     in the small cases for t=4
ColbournKeriStarters <- list(
  '4'=c(CS_CKStarters3[c('12','20')],
        CS_CKStarters1['21'],
        CS_CKStarters3['22'],
        CS_CKStarters2['24'], ## also for 23
        CS_CKStarters3[c('25','26','27')], ## have two constant rows
        CS_CKStarters2[-1])
  )
method4 <- paste0("CS_CK", c(3,3,1,3,2,3,3,3,2,2,2,2,2,2,2))
length(ColbournKeriStarters[['4']]); length(method4)

k4 <- as.numeric(names(ColbournKeriStarters[['4']]))
N4 <- sapply(1:length(method4), function(obj) switch(EXPR=method4[obj],
                                                     NA, CS_CK1=2*k4[obj],
                                                     CS_CK2=2*(k4[obj]-1),
                                                     CS_CK3=2*k4[obj]))

## the part based on starter vectors
ColbournKericat4 <- data.frame(
  t=4, k=as.numeric(names(ColbournKeriStarters[['4']])), v=2, N=N4, method=method4,
  starter=paste0("ColbournKeriStarters[['4']][['",k4,"']]"),
  code=paste0(method4, "(", k4, ", starter=ColbournKeriStarters[['4']][['",k4,"']])"))
## later append second part based on deriving strength 5 CAs


####################################### strength 5 #########################
## Lemma 4
## strength 5
## cycvec for q=67, 71, 79, 83 with CS_CK3 - 67, 71, 83 optimal, 79 almost optimal (+1)
## cycvec for q=103 with CS_CK1            - optimal
## cycvec for q=89,97,101 and 107 <= q <= 773 with CS_CK2    - many optimal,
      # two slightly improved by "Paley type (Colburn)" (+2, q=127, 131),
      # various massively improved by "Paley type (Colburn)" (almost halved)
            ## (q=379 to q=421 and 433 to 491)
      # various massively improved by "Cyclotomy (Colbourn)"
            ## (q=431 (Cyclotomy q=433, construction 2),
            ## q=521 to 607: cyclotomy construction 1, same q for up to q columns,
            ##      next q for k=q+1 columns
      # various slightly or more strongly improved by SBSTT (TJ-AG) (+1 (q=113),
            ## +3 (q=137), +5 (q=139), +9 (q=149),
            ## +10 (q=151, 157), +16 (q=163, 167, 181), +19 (q=173), +26 (q=179),
            ## +28 (q=181), +36 (q=191), +29 (q=193)
      # various massively improved by derive from strength 6
            ## from q=613 (for k=614, for k=613, still cyclotomy)
## derived from strength 6
## k=q for q in 359, 431, 463, 467, 487, 491, 499, 503, 509, 521,
## k=q-1 for q=379

# ## inspect the strength 5 constructions for q=107 to 773
# x <- as.numeric(strsplit("107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281 283 293 307 311 313 317 331 337 347 349 353 359 367
#   373 379 383 389 397 401 409 419 421 431 433 439 443 449 457 461 463 467 479 487 491 499 503
#   509 521 523 541 547 557 563 569 571 577 587 593 599 601 607 613 617 619 631 641 643 647 653
#   659 661 673 677 683 691 701 709 719 727 733 739 743 751 757 761 769 773", " ", fixed=TRUE)[[1]])
# x <- x[!is.na(x)]
# eCANvec <- Vectorize(eCAN, vectorize.args="k", SIMPLIFY = FALSE)
# cbind(k=x+1, N=2*x, do.call(rbind, t(eCANvec(5, x+1, 2))), do.call(rbind, t(eCANvec(5, x, 2))))

## cycvec for q=67, 71, 79, 83 with CS_CK3 - 67, 71, 83 optimal, 79 almost optimal (+1)
## cycvec for q=103 with CS_CK1            - optimal
## cycvec for q=89,97,101 and 107 <= q <= 773 with CS_CK2    - many optimal,

method5 <- paste0("CS_CK", c(3,3,3,3,2,2,2,1,rep(2, 110)))
q5 <- c(67, 71, 79, 83, 89,97, 101, 103,
        primes::primes[primes::primes>=107 & primes::primes<=773])
k5 <- sapply(1:length(method5), function(obj) switch(EXPR=method5[obj],
                                                     NA, CS_CK1=q5[obj],
                                                     CS_CK2=q5[obj]+1,
                                                     CS_CK3=q5[obj]+1))
N5 <- sapply(1:length(method5), function(obj) switch(EXPR=method5[obj],
                                                     NA, CS_CK1=2*k5[obj],
                                                     CS_CK2=2*(k5[obj]-1),
                                                     CS_CK3=2*k5[obj]))

## the part based on starter vectors
ColbournKericat5 <- data.frame(
  t=5, k=k5, v=2, N=N5, method=method5,
  starter=paste0("cycvec(",q5,")"),
  code=paste0(method5, "(k=", k5, ", starter=cycvec(2,", q5,"))"))
## later append second part based on deriving strength 6 CAs

#################################### strength 4 revisited #################################
## now append to ColbournKericat4 the second part based on deriving strength 5 CAs
ColbournKericat4 <- rbind(ColbournKericat4,
                          data.frame(t=4, k=ColbournKericat5$k-1,
                                     v=2, N=ColbournKericat5$N%/%2, ## all original N are even
                                     method="derive",
                                     starter=ColbournKericat5$starter,
                                     code=paste0("derive(",
                                                        ColbournKericat5$code, ")"))
                          )

#################################### strength 6 ###########################################
## Lemma 5
## strength 6
## cycvec for q= 359, 503 with CS_CK3 - optimal
## cycvec for q= 379 with CS_CK1      - optimal
## cycvec for q= 431, 463, 467, 487, 491, 499, 509, 521 with CS_CK2 - optimal

## Lemma 6: derived from one higher strength
## strength 4 for k=q columns for primes from 67 to 773, except for 73 and 103
## strength 4 for k=102 columns for q=103
## strength 5 for 359, 431, 463, 467, 487, 491, 499, 503, 509, 521 columns
## strength 5 for 378 columns from prime 379 (method1)

method6 <- paste0("CS_CK", c(3, 1, 2, 2, 2, 2, 2, 2, 3, 2, 2))
q6 <- c(359, 379, 431, 463, 467, 487, 491, 499, 503, 509, 521)
k6 <- sapply(1:length(method6), function(obj) switch(EXPR=method6[obj],
                                                     NA, CS_CK1=q6[obj],
                                                     CS_CK2=q6[obj]+1,
                                                     CS_CK3=q6[obj]+1))
N6 <- sapply(1:length(method6), function(obj) switch(EXPR=method6[obj],
                                                     NA, CS_CK1=2*k6[obj],
                                                     CS_CK2=2*(k6[obj]-1),
                                                     CS_CK3=2*k6[obj]))

ColbournKericat6 <- data.frame(
  t=6, k=k6, v=2, N=N6, method=method6,
  starter=paste0("cycvec(2,",q6,")"),
  code=paste0(method6, "(k=", k6, ", starter=cycvec(2,", q6,"))")
)
#################################### strength 5 revisited #################################
## now append to ColbournKericat5 the second part based on deriving strength 6 CAs
ColbournKericat5 <- rbind(ColbournKericat5,
                          data.frame(t=5, k=ColbournKericat6$k-1,
                                     v=2, N=ColbournKericat6$N%/%2,
                                     method="derive",
                                     starter=ColbournKericat6$starter,
                                     code=paste0("derive(",
                                                 ColbournKericat6$code, ")"))
)

ColbournKeriCombis <- rbind(ColbournKericat4, ColbournKericat5, ColbournKericat6)
eCANvec <- Vectorize(eCAN, vectorize.args=c("t","k"), SIMPLIFY = FALSE)
eCAKvec <- Vectorize(eCAK, vectorize.args=c("t","N"), SIMPLIFY = FALSE)
compare <- cbind(ColbournKeriCombis[,1:5], do.call(rbind, t(eCANvec(ColbournKeriCombis$t,
                                                         ColbournKeriCombis$k, 2))),
                             do.call(rbind, t(eCAKvec(ColbournKeriCombis$t,
                                                      ColbournKeriCombis$N, 2))))
diff <- compare$N-compare$CAN
fivenum(diff)
ratio <- compare$N/compare$CAN

## verify the dimensions of the construction
for (i in 1:nrow(ColbournKeriCombis)){
  ## 273 rows
  print(i)
  hilf <- dim(eval(parse(text=ColbournKeriCombis$code[i])))
  stopifnot(all(hilf==c(ColbournKeriCombis$N[i], ColbournKeriCombis$k[i])))
}

## verify the strengths for constructions with choose(k,t)<=10000000
## and rows 134 and 135dimensions of the construction
for (i in c(1:25, 134)){
  print(i)
  hilf <- coverage(eval(parse(text=ColbournKeriCombis$code[i])), t=ColbournKeriCombis$t[i])
  stopifnot(all(hilf==1))
}

## might add sample verifications


