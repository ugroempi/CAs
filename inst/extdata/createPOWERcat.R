## the data.frame powerCTcat can be used for implementing the construction
## since version 0.19, it is limited to cases with at most 1 Mio runs
## this removes the need for iterative updating, because there are no powerCT
##      construction ingredients any more

## remove N-CT and "Arc" and "S"
whichCTnormal <- setdiff(grep("CT", colbournBigFrame$Source),
                         c(grep("N-CT", colbournBigFrame$Source),
                           grep("Arc", colbournBigFrame$Source),
                           grep("S", colbournBigFrame$Source)))
powerCTcat <- colbournBigFrame[whichCTnormal,]
powerCTcat <- powerCTcat[which(powerCTcat$N<=1000000),]
N <- powerCTcat$N
k <- powerCTcat$k
v <- powerCTcat$v
plus <- sapply(strsplit(powerCTcat$Source,"+", fixed=TRUE), function(obj) ifelse(length(obj)==1, 0, obj[2]))
plus <- as.numeric(gsub("c","",plus))
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

## need to use a bestN that prioritizes constant rows
## past changes to the global function bestN have deteriorated run sizes for the powerCT construction
## CS_CK instead of PALEY in constr1 and constr2 harmful - prevent by moving it behind PALEY
## Paley instead of CK_doublingCA harmful, prevent by moving CK_doublingCA before PALEY
## fuseBushtCA instead of SCA_Busht harmful, as it does not know that there are q constant rows
##            (though practically irrelevant for the large cases it affects, because these are not doable,
##             at least not on my machine)
## NIST instead of DWYER in constr2 and constr3 harmful

#require(CAs)  # version 0.18 with naively updated powerCTcat

powerCTcat$constr4 <- powerCTcat$constr3 <- powerCTcat$constr2 <- powerCTcat$constr1 <- ""
powerCTcat$c4 <- powerCTcat$c3 <- powerCTcat$c2 <- powerCTcat$c1 <- 1  ## minimum achievable, increase via knowledge
powerCTcat$N4 <- powerCTcat$N3 <- powerCTcat$N2 <- powerCTcat$N1 <- NA
powerCTcat$chi <- NA
powerCTcat$N1 <- mapply(bestN, powerCTcat$t, powerCTcat$w1, powerCTcat$v)
powerCTcat$constr1 <- mapply(function(t,k,v) names(bestN(t,k,v)), powerCTcat$t, powerCTcat$w1, powerCTcat$v)
powerCTcat <- powerCTcat[which(!sapply(powerCTcat$constr1, is.null)),] ## exclude rows for which first array has no implemented construction
powerCTcat$constr1 <- unlist(powerCTcat$constr1)
pick <- which(!is.na(powerCTcat$w2))
powerCTcat$N2[pick] <- mapply(bestN, powerCTcat$t[pick], powerCTcat$w2[pick], powerCTcat$v[pick])
powerCTcat$constr2[pick] <- mapply(function(t,k,v) names(bestN(t,k,v)),
                                   powerCTcat$t[pick], powerCTcat$w2[pick], powerCTcat$v[pick])
powerCTcat$constr2 <- unlist(powerCTcat$constr2)
pick <- which(!is.na(powerCTcat$w3))
powerCTcat$N3[pick] <- mapply(bestN, powerCTcat$t[pick], powerCTcat$w3[pick], powerCTcat$v[pick])
powerCTcat$constr3[pick] <- mapply(function(t,k,v) names(bestN(t,k,v)),
                                   powerCTcat$t[pick], powerCTcat$w3[pick], powerCTcat$v[pick])
powerCTcat$constr3 <- unlist(powerCTcat$constr3)
pick <- which(!is.na(powerCTcat$w4))
powerCTcat$N4[pick] <- mapply(bestN, powerCTcat$t[pick], powerCTcat$w4[pick], powerCTcat$v[pick])
powerCTcat$constr4[pick] <- mapply(function(t,k,v) names(bestN(t,k,v)),
                                   powerCTcat$t[pick], powerCTcat$w4[pick], powerCTcat$v[pick])
powerCTcat$constr4 <- unlist(powerCTcat$constr4)

## replace Paley with miscCA for constr2 and constr3 where w2 and w3 are 8, respectively
## two constant rows instead of one
## modified August 15, 2025
powerCTcat$constr2[which(powerCTcat$w2==8 & powerCTcat$v==2)] <- "miscCA"
powerCTcat$c2[which(powerCTcat$w2==8 & powerCTcat$v==2)] <- 2

powerCTcat$constr3[which(powerCTcat$w3==8 & powerCTcat$v==2)] <- "miscCA"
powerCTcat$c3[which(powerCTcat$w3==8 & powerCTcat$v==2)] <- 2

## information on constant rows
## SCA_Busht has v constant rows, as always v columns only
powerCTcat$c1[which((powerCTcat$constr1=="SCA_Busht" | powerCTcat$constr1=="fuseBushtCA") & powerCTcat$w1<=powerCTcat$v)] <-
  powerCTcat$v[which((powerCTcat$constr1=="SCA_Busht" | powerCTcat$constr1=="fuseBushtCA") & powerCTcat$w1<=powerCTcat$v)]
powerCTcat$c2[which((powerCTcat$constr2=="SCA_Busht" | powerCTcat$constr2=="fuseBushtCA") & powerCTcat$w2<=powerCTcat$v)] <-
  powerCTcat$v[which((powerCTcat$constr2=="SCA_Busht" | powerCTcat$constr2=="fuseBushtCA") & powerCTcat$w2<=powerCTcat$v)]
powerCTcat$c3[which((powerCTcat$constr3=="SCA_Busht" | powerCTcat$constr3=="fuseBushtCA") & powerCTcat$w3<=powerCTcat$v)] <-
  powerCTcat$v[which((powerCTcat$constr3=="SCA_Busht" | powerCTcat$constr3=="fuseBushtCA") & powerCTcat$w3<=powerCTcat$v)]
powerCTcat$c4[which((powerCTcat$constr4=="SCA_Busht" | powerCTcat$constr4=="fuseBushtCA") & powerCTcat$w4<=powerCTcat$v)] <-
  powerCTcat$v[which((powerCTcat$constr4=="SCA_Busht" | powerCTcat$constr4=="fuseBushtCA") & powerCTcat$w4<=powerCTcat$v)]
## t=3: TJ 15 runs (constr1) has 2 constant rows, TJ 10 runs (constr2) as well
## t=6: TJ has only 1 constant run in all cases
powerCTcat$c1[which(powerCTcat$constr1=="TJ" & powerCTcat$t == 3)] <- 2
powerCTcat$c2[which(powerCTcat$constr2=="TJ" & powerCTcat$t == 3)] <- 2
powerCTcat$c4[which((powerCTcat$constr4=="SCA_Busht" | powerCTcat$constr4=="fuseBushtCA") & powerCTcat$w4<=powerCTcat$v)] <-
  powerCTcat$v[which((powerCTcat$constr4=="SCA_Busht" | powerCTcat$constr4=="fuseBushtCA") & powerCTcat$w4<=powerCTcat$v)]
powerCTcat$c2[which(powerCTcat$constr2=="FullFactorial")] <- powerCTcat$v[which(powerCTcat$constr2=="FullFactorial")]
powerCTcat$c2[which(powerCTcat$constr2=="CK_doublingCA")] <- powerCTcat$v[which(powerCTcat$constr2=="CK_doublingCA")]
powerCTcat$c3[which(powerCTcat$constr3=="CK_doublingCA")] <- powerCTcat$v[which(powerCTcat$constr3=="CK_doublingCA")]

## NIST is not touched, as the arrays are large
##    the arrays do not have systematic constant rows,
##    i.e., they must be treated with one_is_enough by maxconstant
## DWYER sometimes has constant rows, dwyerCA moves them to the top
##    check if they are there - and modify cj-values

## PALEY and CS_CK have to be checked for each case
for (i in 1:nrow(powerCTcat)){
  print(i)
  constr1 <- powerCTcat$constr1[i]
  constr2 <- powerCTcat$constr2[i]
  constr3 <- powerCTcat$constr3[i]
  constr4 <- powerCTcat$constr4[i]
  t <- powerCTcat$t[i]; v <- powerCTcat$v[i]; k <- powerCTcat$w1[i]
  if (constr1=="DWYER"){
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
  if (constr1 == "CS_CK"){
    k <- powerCTcat$w1[i]
    nc <- ColbournKeriCombis[which(ColbournKeriCombis$t>=t &
        ColbournKeriCombis$k>=k),]$nconst[which.min(ColbournKeriCombis[which(ColbournKeriCombis$t>=t &
                                                                               ColbournKeriCombis$k>=k),]$N)]
    if (nc > 1) powerCTcat$c1[i] <- nc
    if (constr2=="CS_CK"){
      k <- powerCTcat$w2[i]
      nc <- ColbournKeriCombis[which(ColbournKeriCombis$t>=t &
                                       ColbournKeriCombis$k>=k),]$nconst[which.min(ColbournKeriCombis[which(ColbournKeriCombis$t>=t &
                                                                                                    ColbournKeriCombis$k>=k),]$N)]
      if (nc > 1) powerCTcat$c2[i] <- nc
      if (constr3=="CS_CK"){
        k <- powerCTcat$w3[i]
        nc <- ColbournKeriCombis[which(ColbournKeriCombis$t>=t &
                               ColbournKeriCombis$k>=k),]$nconst[which.min(ColbournKeriCombis[which(ColbournKeriCombis$t>=t &
                                                                                  ColbournKeriCombis$k>=k),]$N)]
        if (nc > 1) powerCTcat$c3[i] <- nc
        ## constr4 is never CS_CK
      }
    }
  }
  if (constr1 == "PALEY"){
    k <- powerCTcat$w1[i]
    nc <- PALEYcat[which(PALEYcat$t>=t &
                           PALEYcat$k>=k),]$nconst[which.min(PALEYcat[which(PALEYcat$t>=t &
                                                                              PALEYcat$k>=k),]$N)]
    if (nc > 1) powerCTcat$c1[i] <- nc
    if (constr2=="PALEY"){
      k <- powerCTcat$w2[i]
      nc <- PALEYcat[which(PALEYcat$t>=t &
                             PALEYcat$k>=k),]$nconst[which.min(PALEYcat[which(PALEYcat$t>=t &
                                                                                PALEYcat$k>=k),]$N)]
      if (nc > 1) powerCTcat$c2[i] <- nc
      if (constr3=="PALEY"){
        k <- powerCTcat$w3[i]
        nc <- PALEYcat[which(PALEYcat$t>=t &
                               PALEYcat$k>=k),]$nconst[which.min(PALEYcat[which(PALEYcat$t>=t &
                                                                                  PALEYcat$k>=k),]$N)]
        if (nc > 1) powerCTcat$c3[i] <- nc
        ## constr4 is never PALEY
      }
    }
  }
}

powerCTcat$chi <- pmax(0, powerCTcat$v -
                  powerCTcat$u1*(powerCTcat$v - ifelse(is.na(powerCTcat$c1), 1, powerCTcat$c1)) -
                  powerCTcat$u2*(powerCTcat$v - ifelse(is.na(powerCTcat$c2), 1, powerCTcat$c2)) -
                  powerCTcat$u3*(powerCTcat$v - ifelse(is.na(powerCTcat$c3), 1, powerCTcat$c3)) -
                  powerCTcat$u4*(powerCTcat$v - ifelse(is.na(powerCTcat$c4), 1, powerCTcat$c4)))

powerCTcat$claimedN <- powerCTcat$N
powerCTcat$claimedk <- powerCTcat$k

## this seems correct for standard cases
## cases with massive deviations from claimed N: "+" in the construction, or poor ingredients
powerCTcat$N <- powerCTcat$chi + (powerCTcat$N1-powerCTcat$c1)*powerCTcat$u1 +
  ifelse(is.na((powerCTcat$N2-powerCTcat$c2)*powerCTcat$u2), 0, (powerCTcat$N2-powerCTcat$c2)*powerCTcat$u2) +
  ifelse(is.na((powerCTcat$N3-powerCTcat$c3)*powerCTcat$u3), 0, (powerCTcat$N3-powerCTcat$c3)*powerCTcat$u3) +
  ifelse(is.na((powerCTcat$N4-powerCTcat$c4)*powerCTcat$u4), 0, (powerCTcat$N4-powerCTcat$c4)*powerCTcat$u4)
## also provide real k
homog <- which(is.na(powerCTcat$w2))
powerCTcat$k[homog] <- powerCTcat$Pbase[homog]^powerCTcat$expon[homog]

table(homog=(1:nrow(powerCTcat)) %in% homog, plus=(1:nrow(powerCTcat)) %in% grep("+", powerCTcat$Source, fixed=TRUE))
## all nonhomogeneous real k equal claimed k
powerCTcat$k[setdiff(1:nrow(powerCTcat), homog)] <- powerCTcat$claimedk[setdiff(1:nrow(powerCTcat), homog)]
table(powerCTcat$k - powerCTcat$claimedk)
## +2167 because of the 10000 cap in the tables
## no longer relevant because of 1000000 cap on claimedN
## - x because of powerCTcat[(1:nrow(powerCTcat)) %in% grep("+", powerCTcat$Source, fixed=TRUE),]
boxplot(N - claimedN ~ constr1, data=powerCTcat, horizontal=TRUE, las=1, ylab="",
        subset=!constr1 %in% c("NIST", "DWYER"))

rownames(powerCTcat) <- NULL

## no power constructions involved
table(powerCTcat$constr1)
table(powerCTcat$constr2)
table(powerCTcat$constr2)
table(powerCTcat$constr4)

which(mapply(N_powerCT, powerCTcat$t, powerCTcat$w1, powerCTcat$v)<powerCTcat$N1)
which(mapply(N_powerCT, powerCTcat$t, powerCTcat$w2, powerCTcat$v)<powerCTcat$N2)
which(mapply(N_powerCT, powerCTcat$t, powerCTcat$w3, powerCTcat$v)<powerCTcat$N3)
which(mapply(N_powerCT, powerCTcat$t, powerCTcat$w4, powerCTcat$v)<powerCTcat$N4)

fivenum(powerCTcat$N)

## save
save(powerCTcat, file="D:/rtests/CAs/data/powerCTcat.rda",
     compress="xz")

load("D:/rtests/CAsbak/powerCTcat_v0.21.rda") ## has powerCT_alt
#load("D:/rtests/CAsbak/powerCTcat_v0.18.rda") ## has v18CT
## the cap on claimedN can be omitted in the future
merged <- merge(v18CT[which(v18CT$claimedN<=1000000),], powerCTcat, by=c("t","k","v","Source"))
merged[merged$N.x<merged$N.y,]
dim(merged[merged$N.x>merged$N.y,])
table(merged[merged$N.x>merged$N.y,]$constr1.x)

Ns(2,20000,6)
Ns(4, 2197,2) ## this uses a strength 5 construction
              ## could improve powerCTcat,
              ## but it is simpler to let N_powerCT use
              ##    all constructions with strength at least t
## sanity checks
system.time(aus <- powerCA(6, 361, 3))  ## smallest N for 3-level
dim(aus)
