## the data.frame powerCTcat can be used for implementing the construction

## remove N-CT and "Arc" and "S"
whichCTnormal <- setdiff(grep("CT", colbournBigFrame$Source),
                         c(grep("N-CT", colbournBigFrame$Source),
                           grep("Arc", colbournBigFrame$Source),
                           grep("S", colbournBigFrame$Source)))
powerCTcat <- colbournBigFrame[whichCTnormal,]
N <- powerCTcat$N
k <- powerCTcat$k
v <- powerCTcat$v
Pbase <- as.numeric(substr(sapply(strsplit(powerCTcat$Source, "^", fixed=TRUE),
                                  function(obj) obj[1]), 9, 99))
exponent <- as.numeric(sapply(strsplit(
              sapply(strsplit(
                sapply(strsplit(
                  sapply(strsplit(
                    powerCTcat$Source, "^", fixed=TRUE), function(obj) obj[2]),
                      "T", fixed=TRUE), function(obj) obj[1]),
                      ",", fixed=TRUE), function(obj) obj[1]),
                      "+", fixed=TRUE), function(obj) obj[1]))
powerCTcat[!k <= Pbase^exponent,] ## only homogeneous DHFs with "+ something"

powerCTcat <- cbind(powerCTcat, Pbase=Pbase, expon=exponent)
powerCTcat$OAforDHF <- ""
powerCTcat[which(!powerCTcat$Pbase %in% primedat$q),]
powerCTcat[which(powerCTcat$Pbase==12 & powerCTcat$expon==3),]$OAforDHF <- "miscCA(3,6,12)"
powerCTcat[which(powerCTcat$Pbase==12 & powerCTcat$expon==2),]$OAforDHF <- "miscCA(2,7,12)"
powerCTcat[which(powerCTcat$Pbase>=68 & !powerCTcat$Pbase %in% primedat$q),]$OAforDHF <- "to be found"
powerCTcat[powerCTcat$Pbase %in% primedat$q,]$OAforDHF <- paste0("SCA_Busht(",
      powerCTcat[powerCTcat$Pbase %in% primedat$q,]$Pbase,
      ",",powerCTcat[powerCTcat$Pbase %in% primedat$q,]$expon,")")
## remove cases that cannot be covered at present
powerCTcat <- powerCTcat[which(!powerCTcat$OAforDHF=="to be found"), ]

powerCTcat$M <- sapply(1:nrow(powerCTcat), function(i){
  print(i)
  oatext <- powerCTcat$OAforDHF[i]
  toa <- powerCTcat$expon[i]
  tdhf <- powerCTcat$t[i]
  vdhf <- powerCTcat$v[i]
  OA <- eval(parse(text=oatext))
  nrow(createDHF(OA, toa, tdhf, vdhf))
})

## check plausibility
TuranVec <- Vectorize(Turan)
Mmin <- (powerCTcat$expon-1)*TuranVec(powerCTcat$t, powerCTcat$v) + 1
any(!Mmin==powerCTcat$M)

## initialize info on the required ingredient OAs
powerCTcat$w1 <- NA; powerCTcat$u1 <- 0
powerCTcat$w2 <- NA; powerCTcat$u2 <- 0
powerCTcat$w3 <- NA; powerCTcat$u3 <- 0
powerCTcat$w4 <- NA; powerCTcat$u4 <- 0

## remove first eight characters (Power CT) from Source entry
hilf <- substr(powerCTcat$Source, 9, 99)
## remove the "c" (whose meaning has not yet been found)
hilf <- gsub("c", "", hilf, fixed=TRUE)
## split at the "T", which indicates numbers of levels to remove
hilf <- strsplit(hilf, "T", fixed=TRUE)

## every array has at least one DHF row in Pbase levels
## treat first DHF row for all cases
powerCTcat$w1 <- powerCTcat$Pbase
powerCTcat$u1 <- powerCTcat$M - (lengths(hilf) - 1)  ## remove number of DHF rows with fewer entries

## treat arrays with *at least* one reduced row
one <- which(lengths(hilf) >= 2)
powerCTcat$w2[one] <- powerCTcat$w1[one] - as.numeric(sapply(hilf[one],
                                                         function(obj) obj[2]))
powerCTcat$u2[one] <- 1 # temporary, increase, if more of the same are found

## at least two reduced rows
## two T
two <- which(lengths(hilf) >= 3)
## same reduction as first case
equal <- two[which(sapply(hilf[two], function(obj) obj[2]==obj[3]))]
## new reduction
other <- setdiff(two, equal)
powerCTcat$w3[other] <- powerCTcat$w1[other] - as.numeric(sapply(hilf[other],
                                                             function(obj) obj[3]))
powerCTcat$u3[other] <- 1 # temporary, increase, if more of the same are found
powerCTcat$u2[equal] <- 2 # w3 and k3 still NA

## at least three  reduced rows (T; more does not occur)
three <- which(lengths(hilf) >= 4)
## same reduction as first row
equal1 <- three[which(sapply(hilf[three], function(obj) obj[2]==obj[4]))]
## not same as first, but same as second
equal2 <- setdiff(three[which(sapply(hilf[three],
                                     function(obj) obj[3]==obj[4]))], equal1)
## new reduction
otherthree <- setdiff(three, c(equal1, equal2, equal))
othertwo <- setdiff(three, c(equal1, equal2, otherthree))

powerCTcat$w4[otherthree] <- powerCTcat$w1[otherthree] -
         as.numeric(sapply(hilf[otherthree], function(obj) obj[4]))
powerCTcat$u4[otherthree] <- 1
powerCTcat$w3[othertwo] <- powerCTcat$w1[othertwo] -
         as.numeric(sapply(hilf[othertwo], function(obj) obj[4]))
powerCTcat$u3[othertwo] <- 1
powerCTcat$u2[equal1] <- 3
powerCTcat$u3[equal2] <- 2   ## not equal to first, but equal to second

## to automate the designs: bestCA(t, wj, v)
## BE ALERT WHETHER THE c in the source has an impact somehow

## determine N1 to N4 initially without powerCT,
## as the N in the Colbourn tables are far too optimistic
## rerun later with realistic N in order to see whether there are improvements

powerCTcat$constr4 <- powerCTcat$constr3 <- powerCTcat$constr2 <- powerCTcat$constr1 <- ""
powerCTcat$c4 <- powerCTcat$c3 <- powerCTcat$c2 <- powerCTcat$c1 <- 1  ## minimum achievable, increase via knowledge
powerCTcat$N4 <- powerCTcat$N3 <- powerCTcat$N2 <- powerCTcat$N1 <- NA
powerCTcat$chi <- NA
powerCTcat$N1 <- mapply(bestN, powerCTcat$t, powerCTcat$w1, powerCTcat$v, MoreArgs = list(exclude="powerCT"))
powerCTcat$constr1 <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat$t, powerCTcat$w1, powerCTcat$v, MoreArgs = list(exclude="powerCT"))
powerCTcat <- powerCTcat[which(!sapply(powerCTcat$constr1, is.null)),] ## exclude rows for which first array has no implemented construction
powerCTcat$constr1 <- unlist(powerCTcat$constr1)
pick <- which(!is.na(powerCTcat$w2))
powerCTcat$N2[pick] <- mapply(bestN, powerCTcat$t[pick], powerCTcat$w2[pick], powerCTcat$v[pick], MoreArgs = list(exclude="powerCT"))
powerCTcat$constr2[pick] <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat$t[pick], powerCTcat$w2[pick], powerCTcat$v[pick], MoreArgs = list(exclude="powerCT"))
powerCTcat$constr2 <- unlist(powerCTcat$constr2)
pick <- which(!is.na(powerCTcat$w3))
powerCTcat$N3[pick] <- mapply(bestN, powerCTcat$t[pick], powerCTcat$w3[pick], powerCTcat$v[pick], MoreArgs = list(exclude="powerCT"))
powerCTcat$constr3[pick] <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat$t[pick], powerCTcat$w3[pick], powerCTcat$v[pick], MoreArgs = list(exclude="powerCT"))
powerCTcat$constr3 <- unlist(powerCTcat$constr3)
pick <- which(!is.na(powerCTcat$w4))
powerCTcat$N4[pick] <- mapply(bestN, powerCTcat$t[pick], powerCTcat$w4[pick], powerCTcat$v[pick], MoreArgs = list(exclude="powerCT"))
powerCTcat$constr4[pick] <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat$t[pick], powerCTcat$w4[pick], powerCTcat$v[pick], MoreArgs = list(exclude="powerCT"))
powerCTcat$constr4 <- unlist(powerCTcat$constr4)

## information on constant rows
## SCA_Busht has v constant rows, as always v columns only
table(unlist(powerCTcat[,c("constr1","constr2","constr3", "constr4")]), rep(1:4, each=116))
powerCTcat$c1[which(powerCTcat$constr1 %in% c("SCA_Busht", "FullFactorial"))] <-
  powerCTcat$v[which(powerCTcat$constr1 %in% c("SCA_Busht", "FullFactorial"))]
powerCTcat$c2[which(powerCTcat$constr2 %in% c("SCA_Busht", "FullFactorial"))] <-
  powerCTcat$v[which(powerCTcat$constr2 %in% c("SCA_Busht", "FullFactorial"))]

## CKRS (single construction) has one constant row only --> keep unchanged
## and tj
## tj 15 runs (constr1) has 2 constant rows, tj12 runs (constr2) 12 doesn't
## for Pbase 19 only 1 constant run

## inspect TJ
powerCTcat[which(powerCTcat$constr1=="TJ" | powerCTcat$constr2=="TJ"),]
## maxconstant(tjCA(3, 5, 2)) ## permits 2 constant rows constr2, N2=10 for t=3, v=2, k=605 and k=6655
## maxconstant(tjCA(3, 12, 2)) ## permits 2 constant rows constr1, N1=12 for t=3, v=2, k=1452 and 1584
powerCTcat$c1[which(powerCTcat$constr1=="TJ" & powerCTcat$Pbase==12)] <- 2
powerCTcat$c2[which(powerCTcat$constr2=="TJ" & powerCTcat$Pbase==11)] <- 2

## inspect CK_doublingCA
powerCTcat[which(powerCTcat$constr2=="CK_doublingCA" | powerCTcat$constr3=="CK_doublingCA"),]
## it's the same in all cases: CK_doublingCA(4,2), two constant rows (the eight run design
##    with two constant rows)
powerCTcat$c2[which(powerCTcat$constr2=="CK_doublingCA")] <- 2
powerCTcat$c3[which(powerCTcat$constr3=="CK_doublingCA")] <- 2

## NIST is not touched, as the arrays are large
##    the arrays do not have systematic constant rows,
##    i.e., they must be treated with one_is_enough by maxconstant
## DWYER sometimes has constant rows, dwyerCA moves them to the top
##    check if they are there - and modify cj-values
## PALEY has to be checked for each case

## treat DWYER
DWYERpick <- which(powerCTcat$constr1=="DWYER" |
                     powerCTcat$constr2=="DWYER" |
                     powerCTcat$constr3=="DWYER")
powerCTcat[DWYERpick,]

for (i in DWYERpick){
  print(i)
  constr1 <- powerCTcat$constr1[i]
  constr2 <- powerCTcat$constr2[i]
  constr3 <- powerCTcat$constr3[i]
  t <- powerCTcat$t[i]; v <- powerCTcat$v[i]
  if (constr1=="DWYER"){
    k <- powerCTcat$w1[i]
    hilf <- dwyerCA(t,k,v)
    nc <- attr(hilf, "nconstant")
    if (!is.null(nc)) powerCTcat$c1[i] <- nc
  }
  if (constr2=="DWYER"){
      k <- powerCTcat$w2[i]
      hilf <- dwyerCA(t,k,v)
      nc <- attr(hilf, "nconstant")
      if (!is.null(nc)) powerCTcat$c2[i] <- nc
  }
  if (constr3=="DWYER"){
        k <- powerCTcat$w3[i]
        hilf <- dwyerCA(t,k,v)
        nc <- attr(hilf, "nconstant")
        if (!is.null(nc)) powerCTcat$c3[i] <- nc
        ## constr4 is never DWYER
  }
}

## treat PALEY
PALEYpick <- which(powerCTcat$constr1=="PALEY" |
                     powerCTcat$constr2=="PALEY" |
                     powerCTcat$constr3=="PALEY")
powerCTcat[PALEYpick,]

for (i in PALEYpick){
  print(i)
  constr1 <- powerCTcat$constr1[i]
  constr2 <- powerCTcat$constr2[i]
  constr3 <- powerCTcat$constr3[i]
  t <- powerCTcat$t[i]; v <- powerCTcat$v[i]

if (constr1 == "PALEY"){
    hilf <- paleyCA(t,k,v)
    nc <- length(attr(maxconstant(hilf, verbose=2), "constant_rows")$row_set_list)
    if (nc > 1) powerCTcat$c1[i] <- nc
}
    if (constr2=="PALEY"){
      k <- powerCTcat$w2[i]
      hilf <- paleyCA(t,k,v)
      nc <- length(attr(maxconstant(hilf, verbose=2), "constant_rows")$row_set_list)
      if (nc > 1) powerCTcat$c2[i] <- nc
      }
      if (constr3=="PALEY"){
        k <- powerCTcat$w3[i]
        hilf <- paleyCA(t,k,v)
        nc <- length(attr(maxconstant(hilf, verbose=2), "constant_rows")$row_set_list)
        if (nc > 1) powerCTcat$c3[i] <- nc
        ## constr4 is never PALEY
      }
}

powerCTcat$chi <- pmax(0, powerCTcat$v -
                  powerCTcat$u1*(powerCTcat$v - ifelse(is.na(powerCTcat$c1), 1, powerCTcat$c1)) -
                  powerCTcat$u2*(powerCTcat$v - ifelse(is.na(powerCTcat$c2), 1, powerCTcat$c2)) -
                  powerCTcat$u3*(powerCTcat$v - ifelse(is.na(powerCTcat$c3), 1, powerCTcat$c3)) -
                  powerCTcat$u4*(powerCTcat$v - ifelse(is.na(powerCTcat$c4), 1, powerCTcat$c4)))

powerCTcat$claimedN <- powerCTcat$N
powerCTcat$claimedk <- powerCTcat$k

## this seems to be still wrong, though I cannot see why
powerCTcat$N <- powerCTcat$chi + (powerCTcat$N1-powerCTcat$c1)*powerCTcat$u1 +
  ifelse(is.na((powerCTcat$N2-powerCTcat$c2)*powerCTcat$u2), 0, (powerCTcat$N2-powerCTcat$c2)*powerCTcat$u2) +
  ifelse(is.na((powerCTcat$N3-powerCTcat$c3)*powerCTcat$u3), 0, (powerCTcat$N3-powerCTcat$c3)*powerCTcat$u3) +
  ifelse(is.na((powerCTcat$N4-powerCTcat$c4)*powerCTcat$u4), 0, (powerCTcat$N4-powerCTcat$c4)*powerCTcat$u4)
## also provide real k
homog <- which(is.na(powerCTcat$w2))
powerCTcat$k[homog] <- powerCTcat$Pbase[homog]^powerCTcat$expon[homog]
## das unten stimmt nicht mehr
table(homog=(1:nrow(powerCTcat)) %in% homog, plus=(1:nrow(powerCTcat)) %in% grep("+", powerCTcat$Source, fixed=TRUE))
## all nonhomogeneous real k equal claimed k
powerCTcat$k[setdiff(1:nrow(powerCTcat), homog)] <- powerCTcat$claimedk[setdiff(1:nrow(powerCTcat), homog)]
table(powerCTcat$k - powerCTcat$claimedk)
## +2167 because of the 10000 cap in the tables
## - x because of powerCTcat[(1:nrow(powerCTcat)) %in% grep("+", powerCTcat$Source, fixed=TRUE),]
boxplot(N - claimedN ~ constr1, data=powerCTcat, horizontal=TRUE, las=1, ylab="", subset=!constr1 %in% c("NIST", "DWYER"))

save(powerCTcat, file="D:/rtests/CAs/data/powerCTcat.rda", compress="xz")

# ## dimensions
# ## did not run completely, but covered each v and each construction variant that was possible with my memory size
# ## could not do the very last final rows (could not even allocate a matrix of NAs of that size)
# for (i in 40:nrow(powerCTcat)){
#   if (i==39) next ## dominated by row 40
#   print(i)
#   print(powerCTcat$constr1[i])
#   hilf <- powerCA(powerCTcat$t[i], powerCTcat$k[i], powerCTcat$v[i])
#   hilf2 <- dim(hilf)
#   if (!(hilf2[1]==powerCTcat$N[i] && (hilf2[2]==powerCTcat$k[i] || powerCTcat$k[i]==10000))) break
#   print("done")
# }
#
#
# ## sanity checks
# ## coverage
# aus <- powerCA(3, 88, 2) ## smallest N and k for 2-level
# coverage(aus, 3)  ## ok
#
# aus <- powerCA(6, 361, 3)  ## smallest N and k for 3-level
# ## this matrix is too large even for a sample check, as R appears to
# ## break down under its sheer size of about 152 MB
# ## (the simple dim statement dim(aus[,pick]) can take quite long at times)
# for (i in 1:100){
#   print(i)
#   pick <- sample(1:361, 20)
#   if (!all(coverage(aus[,pick], 6)==1)) break
# }

###################################################################################################
## now rerun including powerCT

powerCTcat2 <- colbournBigFrame[whichCTnormal,]
N <- powerCTcat2$N
k <- powerCTcat2$k
v <- powerCTcat2$v
Pbase <- as.numeric(substr(sapply(strsplit(powerCTcat2$Source, "^", fixed=TRUE),
                                  function(obj) obj[1]), 9, 99))
exponent <- as.numeric(sapply(strsplit(
  sapply(strsplit(
    sapply(strsplit(
      sapply(strsplit(
        powerCTcat2$Source, "^", fixed=TRUE), function(obj) obj[2]),
      "T", fixed=TRUE), function(obj) obj[1]),
    ",", fixed=TRUE), function(obj) obj[1]),
  "+", fixed=TRUE), function(obj) obj[1]))
powerCTcat2[!k <= Pbase^exponent,] ## only homogeneous DHFs with "+ something"

powerCTcat2 <- cbind(powerCTcat2, Pbase=Pbase, expon=exponent)
powerCTcat2$OAforDHF <- ""
powerCTcat2[which(!powerCTcat2$Pbase %in% primedat$q),]
powerCTcat2[which(powerCTcat2$Pbase==12 & powerCTcat2$expon==3),]$OAforDHF <- "miscCA(3,6,12)"
powerCTcat2[which(powerCTcat2$Pbase==12 & powerCTcat2$expon==2),]$OAforDHF <- "miscCA(2,7,12)"
powerCTcat2[which(powerCTcat2$Pbase>=68 & !powerCTcat2$Pbase %in% primedat$q),]$OAforDHF <- "to be found"
powerCTcat2[powerCTcat2$Pbase %in% primedat$q,]$OAforDHF <- paste0("SCA_Busht(",
                                                                   powerCTcat2[powerCTcat2$Pbase %in% primedat$q,]$Pbase,
                                                                   ",",powerCTcat2[powerCTcat2$Pbase %in% primedat$q,]$expon,")")
## remove cases that cannot be covered at present
powerCTcat2 <- powerCTcat2[which(!powerCTcat2$OAforDHF=="to be found"), ]

powerCTcat2$M <- sapply(1:nrow(powerCTcat2), function(i){
  print(i)
  oatext <- powerCTcat2$OAforDHF[i]
  toa <- powerCTcat2$expon[i]
  tdhf <- powerCTcat2$t[i]
  vdhf <- powerCTcat2$v[i]
  OA <- eval(parse(text=oatext))
  nrow(createDHF(OA, toa, tdhf, vdhf))
})

## check plausibility
TuranVec <- Vectorize(Turan)
Mmin <- (powerCTcat2$expon-1)*TuranVec(powerCTcat2$t, powerCTcat2$v) + 1
any(!Mmin==powerCTcat2$M)

## initialize info on the required ingredient OAs
powerCTcat2$w1 <- NA; powerCTcat2$u1 <- 0
powerCTcat2$w2 <- NA; powerCTcat2$u2 <- 0
powerCTcat2$w3 <- NA; powerCTcat2$u3 <- 0
powerCTcat2$w4 <- NA; powerCTcat2$u4 <- 0

## remove first eight characters (Power CT) from Source entry
hilf <- substr(powerCTcat2$Source, 9, 99)
## remove the "c" (whose meaning has not yet been found)
hilf <- gsub("c", "", hilf, fixed=TRUE)
## split at the "T", which indicates numbers of levels to remove
hilf <- strsplit(hilf, "T", fixed=TRUE)

## every array has at least one DHF row in Pbase levels
## treat first DHF row for all cases
powerCTcat2$w1 <- powerCTcat2$Pbase
powerCTcat2$u1 <- powerCTcat2$M - (lengths(hilf) - 1)  ## remove number of DHF rows with fewer entries

## treat arrays with *at least* one reduced row
one <- which(lengths(hilf) >= 2)
powerCTcat2$w2[one] <- powerCTcat2$w1[one] - as.numeric(sapply(hilf[one],
                                                               function(obj) obj[2]))
powerCTcat2$u2[one] <- 1 # temporary, increase, if more of the same are found

## at least two reduced rows
## two T
two <- which(lengths(hilf) >= 3)
## same reduction as first case
equal <- two[which(sapply(hilf[two], function(obj) obj[2]==obj[3]))]
## new reduction
other <- setdiff(two, equal)
powerCTcat2$w3[other] <- powerCTcat2$w1[other] - as.numeric(sapply(hilf[other],
                                                                   function(obj) obj[3]))
powerCTcat2$u3[other] <- 1 # temporary, increase, if more of the same are found
powerCTcat2$u2[equal] <- 2 # w3 and k3 still NA

## at least three  reduced rows (T; more does not occur)
three <- which(lengths(hilf) >= 4)
## same reduction as first row
equal1 <- three[which(sapply(hilf[three], function(obj) obj[2]==obj[4]))]
## not same as first, but same as second
equal2 <- setdiff(three[which(sapply(hilf[three],
                                     function(obj) obj[3]==obj[4]))], equal1)
## new reduction
otherthree <- setdiff(three, c(equal1, equal2, equal))
othertwo <- setdiff(three, c(equal1, equal2, otherthree))

powerCTcat2$w4[otherthree] <- powerCTcat2$w1[otherthree] -
  as.numeric(sapply(hilf[otherthree], function(obj) obj[4]))
powerCTcat2$u4[otherthree] <- 1
powerCTcat2$w3[othertwo] <- powerCTcat2$w1[othertwo] -
  as.numeric(sapply(hilf[othertwo], function(obj) obj[4]))
powerCTcat2$u3[othertwo] <- 1
powerCTcat2$u2[equal1] <- 3
powerCTcat2$u3[equal2] <- 2   ## not equal to first, but equal to second

## to automate the designs: bestCA(t, wj, v)
## BE ALERT WHETHER THE c in the source has an impact somehow

## determine N1 to N4 initially without powerCT,
## as the N in the Colbourn tables are far too optimistic
## rerun later with realistic N in order to see whether there are improvements

powerCTcat2$constr4 <- powerCTcat2$constr3 <- powerCTcat2$constr2 <- powerCTcat2$constr1 <- ""
powerCTcat2$c4 <- powerCTcat2$c3 <- powerCTcat2$c2 <- powerCTcat2$c1 <- 1  ## minimum achievable, increase via knowledge
powerCTcat2$N4 <- powerCTcat2$N3 <- powerCTcat2$N2 <- powerCTcat2$N1 <- NA
powerCTcat2$chi <- NA
powerCTcat2$N1 <- mapply(bestN, powerCTcat2$t, powerCTcat2$w1, powerCTcat2$v)
powerCTcat2$constr1 <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat2$t, powerCTcat2$w1, powerCTcat2$v)
powerCTcat2 <- powerCTcat2[which(!sapply(powerCTcat2$constr1, is.null)),] ## exclude rows for which first array has no implemented construction
powerCTcat2$constr1 <- unlist(powerCTcat2$constr1)
pick <- which(!is.na(powerCTcat2$w2))
powerCTcat2$N2[pick] <- mapply(bestN, powerCTcat2$t[pick], powerCTcat2$w2[pick], powerCTcat2$v[pick])
powerCTcat2$constr2[pick] <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat2$t[pick], powerCTcat2$w2[pick], powerCTcat2$v[pick])
powerCTcat2$constr2 <- unlist(powerCTcat2$constr2)
pick <- which(!is.na(powerCTcat2$w3))
powerCTcat2$N3[pick] <- mapply(bestN, powerCTcat2$t[pick], powerCTcat2$w3[pick], powerCTcat2$v[pick])
powerCTcat2$constr3[pick] <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat2$t[pick], powerCTcat2$w3[pick], powerCTcat2$v[pick])
powerCTcat2$constr3 <- unlist(powerCTcat2$constr3)
pick <- which(!is.na(powerCTcat2$w4))
powerCTcat2$N4[pick] <- mapply(bestN, powerCTcat2$t[pick], powerCTcat2$w4[pick], powerCTcat2$v[pick])
powerCTcat2$constr4[pick] <- mapply(function(t,k,v, ...) names(bestN(t,k,v, ...)), powerCTcat2$t[pick], powerCTcat2$w4[pick], powerCTcat2$v[pick])
powerCTcat2$constr4 <- unlist(powerCTcat2$constr4)

## information on constant rows
## SCA_Busht has v constant rows, as always v columns only
table(unlist(powerCTcat2[,c("constr1","constr2","constr3", "constr4")]), rep(1:4, each=157))
## this adds 41 CAs that could not be constructed otherwise, with at least 100000 runs in eCAN,
## and MUCH more in real --> not attractive, omit for the moment

## for SCA_Busht in constr2, c2=1 (q+1 columns)
