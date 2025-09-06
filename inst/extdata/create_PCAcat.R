## PCA constructions involving CS constructions with k1=k-1
## and some miscCA arrays
## stated numbers in the comments have not been kept up-to-date
## stages 3 to 5 could be done in a loop

################ stage 1 #############
colnames(CMMSSYCombis)
colnames(LCDSTCombis)
colnames(MeagherStevensCombis)
colnames(miscCAcat)

PCA_cands <- rbind(
  cbind(CMMSSYCombis[,c("v","k","N")], type="CMMSSY"),
  cbind(LCDSTCombis[which(LCDSTCombis$f==1),c("v","k","N")], type="LCDST"),
  cbind(as.data.frame(MeagherStevensCombis[,c("v","k","N")]), type="MS"),
  cbind(miscCAcat[miscCAcat$t==2 & miscCAcat$PCAstatus > 1, c("v", "k", "N")], type="misc")
)

PCA_cands$k1 <- NA
PCA_cands$k1[PCA_cands$type %in% c("CMMSSY", "LCDST", "MS")] <-
  PCA_cands$k[PCA_cands$type %in% c("CMMSSY", "LCDST", "MS")] - 1
PCA_cands$k1[PCA_cands$type=="misc"] <- miscCAcat[miscCAcat$t==2 & miscCAcat$PCAstatus > 1, "PCAstatus"]
PCA_cands$code <- sapply(1:nrow(PCA_cands),
                         function(i) switch(PCA_cands$type[i],
                         "CMMSSY"=paste0("CS_CMMSSY(",PCA_cands$k[i],", ", PCA_cands$v[i], ")"),
                         "LCDST"=paste0("SCA_LCDST(",PCA_cands$k[i],", ", PCA_cands$v[i], ")"),
                         "MS"=paste0("SCA_MS(",PCA_cands$k[i],", ", PCA_cands$v[i], ")"),
                         "misc"=paste0("miscCA(2, ",PCA_cands$k[i],", ", PCA_cands$v[i], ")")))
PCA_cands <- PCA_cands[ord(PCA_cands[,c("v","k")]),]
PCA_cands[PCA_cands$k==7 & PCA_cands$v==5 & PCA_cands$type=="misc",]$code <- "CA_to_PCA(miscCA(2,7,5), tryhard = TRUE)"
PCA_cands[PCA_cands$k==3 & PCA_cands$v==6,]$code <- "CA_to_PCA(miscCA(2,3,6))"
PCA_cands[PCA_cands$k==4 & PCA_cands$v==10,]$code <- "CA_to_PCA(miscCA(2, 4, 10),tryhard = TRUE)"
PCA_cands[PCA_cands$k==7 & PCA_cands$v==12,]$code <- "CA_to_PCA(miscCA(2,7,12))"

###### remove dominated or equally good rows (where N not smaller and neither k nor k1 larger) ######
uniques <- unique(PCA_cands[,c("v","N")])
dim(uniques)
for (i in 1:nrow(uniques)){
  pick <- which(  PCA_cands$v==uniques$v[i] &
                    PCA_cands$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCA_cands$k[pick])
  keepi1 <- which.max(PCA_cands$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCA_cands <- PCA_cands[-remove,]
}
dim(PCA_cands) ## 131 rows from stage1
PCA_cands$stage <- 1

######## remove dominated stage 1 rows, if any ####
remove <- numeric(0)
for (i in 2:nrow(PCA_cands)){
  ki <- PCA_cands$k[i]
  ki1 <- PCA_cands$k1[i]

  vi <- PCA_cands$v[i]
  Ni <- PCA_cands$N[i]
  betteri <- which((PCA_cands[1:(i-1),]$k > ki |
                      PCA_cands[1:(i-1),]$k1 > ki1 |
                     PCA_cands[1:(i-1),]$N < Ni) &
                     PCA_cands[1:(i-1),]$k >= ki &
                     PCA_cands[1:(i-1),]$k1 >= ki1 &
                     PCA_cands[1:(i-1),]$v == vi &
                     PCA_cands[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
if (length(remove)>0)
   PCA_cands <- PCA_cands[-remove,]
dim(PCA_cands) ## 129x7
rownames(PCA_cands) <- NULL

## implement PCA constructions for v >= 3
## with ingredients from CS_CMMSSY, CS_MS, CS_LCDST for f=1, miscCAcat

PCAcat <- PCA_cands
PCAcat <- PCAcat[-(1:nrow(PCAcat)),]

######### stage 2 ####################
## add pairwise combinations among arrays with the same v
for (v in unique(PCA_cands$v)){
  hilf <- PCA_cands[which(PCA_cands$v==v),]
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     k1=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf)){
      k <- hilf[i,"k"]
      k1 <- hilf$k1[i]
      k2 <- k - k1
      typek <- hilf[i,"type"]
      codek <- hilf[i,"code"]
    for (j in i:nrow(hilf)){
      l <- hilf[j,"k"]
      N <- hilf$N[i] + hilf$N[j] - v
      l1 <- hilf$k1[j]
      l2 <- l - l1
      typel <- hilf[j,"type"]
      codel <- hilf[j,"code"]
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N,
                                     type=paste0(typek,"x",typel),
                                       k1=k1neu,
                                     code=paste0("productPCA(",codek, ", ",
                                     codel, ")"),
                                     stage=2L))
    }
  }
  PCAcat <- rbind(PCAcat, addv)
}

rownames(PCAcat) <- NULL
dim(PCAcat)  ## 551 x 7

###### remove duplicated stage 2 cases #####
##  by keeping the one(s) with largest k and/or k1
## and the same or smaller N
uniques <- unique(PCAcat[,c("v","N")])
dim(uniques)  ## 268x2
for (i in 1:nrow(uniques)){
  pick <- which(  PCAcat$v==uniques$v[i] &
                    PCAcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCAcat$k[pick])
  keepi1 <- which.max(PCAcat$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCAcat <- PCAcat[-remove,]
}
dim(PCAcat) ## 275 rows
rownames(PCAcat) <- NULL

######## remove dominated stage 2 rows, if any ####
remove <- numeric(0)
for (i in 1:(nrow(PCAcat)-1)){
  ki <- PCAcat$k[i]
  ki1 <- PCAcat$k1[i]

  vi <- PCAcat$v[i]
  Ni <- PCAcat$N[i]
  betteri <- which((PCAcat[1:(i-1),]$k > ki |
                      PCAcat[1:(i-1),]$k1 > ki1 |
                      PCAcat[1:(i-1),]$N < Ni) &
                     PCAcat[1:(i-1),]$k >= ki &
                     PCAcat[1:(i-1),]$k1 >= ki1 &
                     PCAcat[1:(i-1),]$v == vi &
                     PCAcat[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
if (length(remove) > 0)
   PCAcat <- PCAcat[-remove,]
dim(PCAcat) ## 362x7
rownames(PCAcat) <- NULL
table(PCAcat$stage)

# ###### check coverage for stage 2 ###############
# for (i in 1:nrow(PCAcat)){
#   ## dimensions checked completely,
#   ## coverage only up to approx 335
#   ## should be OK for the others as well
#   ## all principles were now involved
#   print(i)
#   plan <- eval(parse(text=PCAcat$code[i]))
#   if (!nrow(plan)==PCAcat$N[i] && ncol(plan)==PCAcat$k[i]){
#     print("wrong dimensions")
#     break
#   }
#   # if (!all(coverage(plan,2)==1)){
#   #   print("coverage violated")
#   #   break
#   # }
# }

##### add SCA_Bose entries to PCA_cands ###############
## recursiveBose takes care of pairs of them, therefore only now
for (v in unique(PCA_cands$v)){
  if (!v %in% primedat$q) next
  PCA_cands <- rbind(PCA_cands,
                     data.frame(v=v, k=v+1, N=v^2,
                                type="Bose",
                                k1=v, code=paste0("SCA_Bose(", v, ")"), stage=2L))
}
dim(PCA_cands)  ## 138x7

###### remove duplicated stage 1 rows ######
##  by keeping the one(s) with largest k and/or k1
uniques <- unique(PCA_cands[,c("v","N")])
dim(uniques)  ## 135
for (i in 1:nrow(uniques)){
  pick <- which(  PCA_cands$v==uniques$v[i] &
                    PCA_cands$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCA_cands$k[pick])
  keepi1 <- which.max(PCA_cands$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCA_cands <- PCA_cands[-remove,]
}
dim(PCA_cands) ## 136 rows from stage1, part 2

######## remove dominated stage 1 rows, if any ####
remove <- numeric(0)
for (i in 2:nrow(PCA_cands)){
  ki <- PCA_cands$k[i]
  ki1 <- PCA_cands$k1[i]

  vi <- PCA_cands$v[i]
  Ni <- PCA_cands$N[i]
  betteri <- which((PCA_cands[1:(i-1),]$k > ki |
                      PCA_cands[1:(i-1),]$k1 > ki1 |
                      PCAcat[1:(i-1),]$N < Ni) &
                     PCA_cands[1:(i-1),]$k >= ki &
                     PCA_cands[1:(i-1),]$k1 >= ki1 &
                     PCA_cands[1:(i-1),]$v == vi &
                     PCA_cands[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
## nothing redundant
rownames(PCA_cands) <- NULL

# ############################ stage 3 ###################################
for (v in unique(PCAcat$v)){
  ## loop over v
  hilf2 <- PCAcat[which(PCAcat$v==v & PCAcat$N<=20000),]  ## stage 2
  hilf1 <- PCA_cands[which(PCA_cands$v==v),]
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0), type=character(0),
                     k1=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    k1 <- hilf1$k1[i]
    k2 <- k - k1
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      N <- hilf1$N[i] + hilf2$N[j] - v
      l1 <- hilf2$k1[j]
      l2 <- l - l1
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1

      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                                     type=paste0(typek,"x",typel),
                                     code=paste0("productPCA(", codek, ", ",
                                                 codel, ")"),
                                     stage=3L))
    }
  }
  PCAcat <- rbind(PCAcat, addv)
}
dim(PCAcat) ## 2417 rows

###### redo the duplication removal #################
uniques <- unique(PCAcat[,c("v","N")])
dim(uniques)  #636x2

for (i in 1:nrow(uniques)){
  pick <- which(  PCAcat$v==uniques$v[i] &
                    PCAcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCAcat$k[pick])
  keepi1 <- which.max(PCAcat$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCAcat <- PCAcat[-remove,]
}
dim(PCAcat) ## 655 rows, 194 from stage 2
table(PCAcat$stage)
rownames(PCAcat) <- NULL

###### check for dominated rows ###############################################
## remove dominated rows, if any
remove <- numeric(0)
removals <- vector(mode="list")
for (i in (nrow(PCAcat):2)){
  ki <- PCAcat$k[i]
  ki1 <- PCAcat$k1[i]

  vi <- PCAcat$v[i]
  Ni <- PCAcat$N[i]
  betteri <- i + which((PCAcat[1:(i-1),]$k > ki |
                          PCAcat[1:(i-1),]$k1 > ki1 |
                          PCAcat[1:(i-1),]$N < Ni) &
                         PCAcat[1:(i-1),]$k >= ki &
                         PCAcat[1:(i-1),]$k1 >= ki1 &
                         PCAcat[1:(i-1),]$v == vi &
                         PCAcat[1:(i-1),]$N <= Ni)
  if (length(betteri)>0){
    print(i)
    print(betteri)
    remove <- c(remove,i)
    removals[[length(removals)+1]] <- list(better=PCAcat$code[betteri],
                                           removed=PCAcat$code[i])
  }
}
if (length(remove) > 0)
  PCAcat <- PCAcat[-remove,]
dim(PCAcat) # 618x7
## limit to k<=20000
## PCAcat <- PCAcat[PCAcat$k<=20000,]
dim(PCAcat) # 547x7
rownames(PCAcat) <- NULL

# ###################### stage 4 #############################
## add fourth factors among arrays with the same v
for (v in unique(PCAcat$v)){
  ## loop over v
  hilf2 <- PCAcat[which(PCAcat$v==v & PCAcat$N<=20000 & PCAcat$stage==3),]
  hilf1 <- PCA_cands[which(PCA_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     k1=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    k1 <- hilf1$k1[i]
    k2 <- k - k1
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      N <- hilf1$N[i] + hilf2$N[j] - v
      l1 <- hilf2$k1[j]
      l2 <- l - l1
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                                     type=paste0(typek,"x",typel),
                                     code=paste0("productPCA(", codek, ", ",
                                                 codel, ")"),
                                     stage=4L))
    }
  }
  PCAcat <- rbind(PCAcat, addv)
}
dim(PCAcat)  ## 3627 rows
table(PCAcat$stage)
## PCAcat <- PCAcat[PCAcat$k<=20000,]
dim(PCAcat)  ## 1675 rows

###### redo the duplication removal within stage 4 ####################
uniques <- unique(PCAcat[,c("v","N")])
dim(uniques)  # 669x2

for (i in 1:nrow(uniques)){
  pick <- which(  PCAcat$v==uniques$v[i] &
                    PCAcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCAcat$k[pick])
  keepi1 <- which.max(PCAcat$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCAcat <- PCAcat[-remove,]
}
dim(PCAcat) ## 703 rows, 189 from stage 2, 257 from stage 3, 257 from stage 4
table(PCAcat$stage)
rownames(PCAcat) <- NULL

###### check for dominated rows ###############################################
## remove dominated rows, if any
remove <- numeric(0)
removals <- vector(mode="list")
for (i in (nrow(PCAcat):2)){
  ki <- PCAcat$k[i]
  ki1 <- PCAcat$k1[i]

  vi <- PCAcat$v[i]
  Ni <- PCAcat$N[i]
  betteri <- i + which((PCAcat[1:(i-1),]$k > ki |
                          PCAcat[1:(i-1),]$k1 > ki1 |
                          PCAcat[1:(i-1),]$N < Ni) &
                         PCAcat[1:(i-1),]$k >= ki &
                         PCAcat[1:(i-1),]$k1 >= ki1 &
                         PCAcat[1:(i-1),]$v == vi &
                         PCAcat[1:(i-1),]$N <= Ni)
  if (length(betteri)>0){
    print(i)
    print(betteri)
    remove <- c(remove,i)
    removals[[length(removals)+1]] <- list(better=PCAcat$code[betteri],
                                           removed=PCAcat$code[i])
  }
}
if (length(remove) > 0)
  PCAcat <- PCAcat[-remove,]
dim(PCAcat) # 618x7

lapply(unique(PCAcat$v), function(obj) fivenum(PCAcat$k[PCAcat$v==obj]))

# ###################### stage 5 #############################
## add fifth factors among arrays with the same v
for (v in unique(PCAcat$v)){
  ## loop over v
  hilf2 <- PCAcat[which(PCAcat$v==v & PCAcat$N<=20000 & PCAcat$stage==4),]
  hilf1 <- PCA_cands[which(PCA_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     k1=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    k1 <- hilf1$k1[i]
    k2 <- k - k1
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      N <- hilf1$N[i] + hilf2$N[j] - v
      l1 <- hilf2$k1[j]
      l2 <- l - l1
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                                     type=paste0(typek,"x",typel),
                                     code=paste0("productPCA(", codek, ", ",
                                                 codel, ")"),
                                     stage=5L))
    }
  }
  PCAcat <- rbind(PCAcat, addv)
}
dim(PCAcat)  ## 2310 rows
table(PCAcat$stage)
## PCAcat <- PCAcat[PCAcat$k<=20000,]
dim(PCAcat)  ## 1311 rows

###### redo the duplication removal within stage 5 ####################
uniques <- unique(PCAcat[,c("v","N")])
dim(uniques)  # 644x2

for (i in 1:nrow(uniques)){
  pick <- which(  PCAcat$v==uniques$v[i] &
                    PCAcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCAcat$k[pick])
  keepi1 <- which.max(PCAcat$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCAcat <- PCAcat[-remove,]
}
dim(PCAcat) ## 686 rows, 122 from stage1, 327 from stages 1 and 2, 549 from 1 to 3, 689 from 1 to 4
table(PCAcat$stage)
rownames(PCAcat) <- NULL

###### check for dominated rows ###############################################
## remove dominated rows, if any
remove <- numeric(0)
removals <- vector(mode="list")
for (i in (nrow(PCAcat):2)){
  ki <- PCAcat$k[i]
  ki1 <- PCAcat$k1[i]

  vi <- PCAcat$v[i]
  Ni <- PCAcat$N[i]
  betteri <- i + which((PCAcat[1:(i-1),]$k > ki |
                          PCAcat[1:(i-1),]$k1 > ki1 |
                          PCAcat[1:(i-1),]$N < Ni) &
                         PCAcat[1:(i-1),]$k >= ki &
                         PCAcat[1:(i-1),]$k1 >= ki1 &
                         PCAcat[1:(i-1),]$v == vi &
                         PCAcat[1:(i-1),]$N <= Ni)
  if (length(betteri)>0){
    print(i)
    print(betteri)
    remove <- c(remove,i)
    removals[[length(removals)+1]] <- list(better=PCAcat$code[betteri],
                                           removed=PCAcat$code[i])
  }
}
if (length(remove) > 0)
  PCAcat <- PCAcat[-remove,]
dim(PCAcat) # 637x7
rownames(PCAcat) <- NULL

# ###################### stage 6 #############################
## add fifth factors among arrays with the same v
for (v in unique(PCAcat$v)){
  ## loop over v
  hilf2 <- PCAcat[which(PCAcat$v==v & PCAcat$N<=20000 & PCAcat$stage==5),]
  hilf1 <- PCA_cands[which(PCA_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     k1=numeric(0), code=character(0),
                     stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    k1 <- hilf1$k1[i]
    k2 <- k - k1
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      N <- hilf1$N[i] + hilf2$N[j] - v
      l1 <- hilf2$k1[j]
      l2 <- l - l1
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                                     type=paste0(typek,"x",typel),
                                     code=paste0("productPCA(", codek, ", ",
                                                 codel, ")"),
                                     stage=6L))
    }
  }
  PCAcat <- rbind(PCAcat, addv)
}
dim(PCAcat)  ## 1520 rows
table(PCAcat$stage)
# at the last stage, keep these, as they may have
#   reasonable sizes
# PCAcat <- PCAcat[PCAcat$k<=20000,]
dim(PCAcat)  ## 822 rows

###### redo the duplication removal within stage 6 ####################
uniques <- unique(PCAcat[,c("v","N")])
dim(uniques)  # 617x2

for (i in 1:nrow(uniques)){
  pick <- which(  PCAcat$v==uniques$v[i] &
                    PCAcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(PCAcat$k[pick])
  keepi1 <- which.max(PCAcat$k1[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    PCAcat <- PCAcat[-remove,]
}
dim(PCAcat) ## 660 rows, 189 from stage2, 218 from stage 3, 116 from 4, 86 from 5
table(PCAcat$stage)
rownames(PCAcat) <- NULL

###### final check for dominated rows ###############################################
###### k1 does not matter any more
## remove dominated rows, if any
remove <- numeric(0)
for (i in 1:nrow(PCAcat)){
  ki <- PCAcat$k[i]
  vi <- PCAcat$v[i]
  Ni <- PCAcat$N[i]
  betteri <- which((PCAcat$k > ki |
                      PCAcat$N < Ni) &
                     PCAcat$v == vi &
                     PCAcat$N <= Ni &
                     PCAcat$k >=ki)

  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
if (length(remove) > 0)
 PCAcat <- PCAcat[-remove,]
dim(PCAcat) # 488x7

PCAcat <- PCAcat[ord(PCAcat[,c("v","k")]),]
rownames(PCAcat) <- NULL

sapply(unique(PCAcat$v), function(obj) fivenum(PCAcat$k[PCAcat$v==obj]))
sapply(unique(PCAcat$v), function(obj) fivenum(PCAcat$N[PCAcat$v==obj]))

## limit so that k=20000 remains covered
remove <- numeric(0)
for (v in unique(PCAcat$v)){
   hilf <- which(PCAcat$v==v & PCAcat$k > 20000)
   remove <- c(remove, setdiff(hilf, hilf[which.min(PCAcat[hilf,]$N)]))
}
PCAcat <- PCAcat[-remove,]
dim(PCAcat)
rownames(PCAcat) <- NULL

table(isbest <- bestNvec(PCAcat$k, PCAcat$v)==PCAcat$N, PCAcat$v)
table(isbest, PCAcat$stage)
table(isopt <- eCANvec(PCAcat$k, PCAcat$v)==PCAcat$N, PCAcat$v)
table(isopt, PCAcat$stage)

######## check dimensions for stages 3 to 5 ################
## stage
table(PCAcat$stage)
for (i in 1:nrow(PCAcat)){
  ## dimensions checked completely,
  print(i)
  plan <- pcaCA(PCAcat$k[i], PCAcat$v[i])
  # plan <- eval(parse(text=PCAcat$code[i]))
  if (!nrow(plan)==PCAcat$N[i] && ncol(plan)==PCAcat$k[i]){
    print("wrong dimensions")
    break
  }
}

# save(PCA_cands, file="D:/rtests/CAs/inst/extdata/PCA_cands.rda", compress="xz")
# save(PCAcat, file="D:/rtests/CAs/data/PCAcat.rda", compress="xz")
