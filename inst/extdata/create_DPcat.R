## CA constructions involving CS constructions with k1=k-1
## and selected arrays from miscCAcat
## stated numbers in the comments have not been kept up-to-date
## stages 3 to 5 could be done in a loop

################ stage 1 #############
colnames(CMMSSYCombis)
colnames(LCDSTCombis)
colnames(MeagherStevensCombis)
colnames(miscCAcat)

mischilfconst <- miscCAcat[which(miscCAcat$PCAstatus<miscCAcat$k &
                                   miscCAcat$nconst > 1 & miscCAcat$t==2),]
mischilfPCA <- miscCAcat[which(miscCAcat$PCAstatus > 1 & miscCAcat$t==2),]
mischilfPCA$k <- mischilfPCA$PCAstatus
mischilfPCA$nconst <- mischilfPCA$v
mischilf <- rbind(mischilfconst, mischilfPCA)
mischilf$PCAstatus<-NULL
mischilf$comment <- NULL
mischilf$hasNA <- NULL
mischilf$fns <- NULL
dim(mischilf)

## remove duplicate
uniques <- unique(mischilf[,c("k","v")])
dim(uniques)
for (i in 1:nrow(uniques)){
  pick <- which(  mischilf$v==uniques$v[i] &
                    mischilf$k==uniques$k[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.min(mischilf$N[pick])
  ## checked manually that this is advantageous for all cases
  ##    also consideting nconst
  ## remove all others
  remove <- pick[-c(keepi)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    mischilf <- mischilf[-remove,]
}
dim(mischilf) ## 18 rows
rownames(mischilf) <- NULL

############### check nconst #######
for (i in 1:nrow(mischilf)){
  print(i)
  nconsti <- mischilf$nconst[i]
  if (!length(unique(maxconstant(miscCA(2,mischilf$k[i],mischilf$v[i]))[nconsti,]))==1) message("problem")
}

## DP_cands analogous to create_PCAcat
##    use the first k1 columns with constant rows
##    bring in mischilf later
DP_cands <- rbind(
  cbind(CMMSSYCombis[,c("v","k","N")], type="CMMSSY"),
  cbind(LCDSTCombis[which(LCDSTCombis$f==1),c("v","k","N")], type="LCDST"),
  cbind(as.data.frame(MeagherStevensCombis[,c("v","k","N")]), type="MS")
)

DP_cands$k <- DP_cands$k - 1 ## reduce to k1
DP_cands$nconst <- DP_cands$v
DP_cands <- rbind(DP_cands, data.frame(v=mischilf$v, k=mischilf$k, N=mischilf$N, type="misc", nconst=mischilf$nconst))

DP_cands$code <- sapply(1:nrow(DP_cands),
                         function(i) switch(DP_cands$type[i],
                         "CMMSSY"=paste0("CS_CMMSSY(",DP_cands$k[i]+1,", ", DP_cands$v[i], ")[,1:", DP_cands$k[i],"]"),
                         "LCDST"=paste0("SCA_LCDST(",DP_cands$k[i]+1,", ", DP_cands$v[i], ")[,1:", DP_cands$k[i],"]"),
                         "MS"=paste0("SCA_MS(",DP_cands$k[i]+1,", ", DP_cands$v[i], ")[,1:", DP_cands$k[i],"]"),
                         "misc"=paste0("miscCA(2, ",DP_cands$k[i],", ", DP_cands$v[i], ")")))
DP_cands <- DP_cands[ord(DP_cands[,c("v","k")]),]
rownames(DP_cands) <- NULL

###### check nconst ##########
## no problems
# for (i in 1:nrow(DP_cands)){
#   print(i)
#   nconsti <- DP_cands$nconst[i]
#   hilf <- maxconstant(eval(parse(text=DP_cands$code[i])))
#   if (!length(unique(hilf[nconsti,]))==1) message("problem")
# }

###### remove dominated or equally good rows (where N not smaller and neither k nor k1 larger) ######
uniques <- unique(DP_cands[,c("v","N")])
dim(uniques)
for (i in 1:nrow(uniques)){
  pick <- which(  DP_cands$v==uniques$v[i] &
                    DP_cands$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DP_cands$k[pick])
  keepi1 <- which.max(DP_cands$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DP_cands <- DP_cands[-remove,]
}
dim(DP_cands) ## 136 rows from stage1
DP_cands$stage <- 1

######## remove dominated stage 1 rows, if any ####
remove <- numeric(0)
for (i in 1:(nrow(DP_cands)-1)){
  ki <- DP_cands$k[i]
  ci <- DP_cands$nconst[i]

  vi <- DP_cands$v[i]
  Ni <- DP_cands$N[i]
  betteri <- which((DP_cands[(i+1):nrow(DP_cands),]$k > ki |
                      DP_cands[(i+1):nrow(DP_cands),]$nconst > ci) &
                     DP_cands[(i+1):nrow(DP_cands),]$k >= ki &
                     DP_cands[(i+1):nrow(DP_cands),]$nconst >= ci &
                     DP_cands[(i+1):nrow(DP_cands),]$v == vi &
                     DP_cands[(i+1):nrow(DP_cands),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
DP_cands <- DP_cands[-remove,]
dim(DP_cands) ## 132x7
rownames(DP_cands) <- NULL

## implement generalized direct product constructions for v >= 3
## with ingredients from CS_CMMSSY, CS_MS, CS_LCDST for f=1, miscCAcat

DPcat <- DP_cands
DPcat <- DPcat[-(1:nrow(DPcat)),]

######### stage 2 ####################
## add pairwise combinations among arrays with the same v
for (v in unique(DP_cands$v)){
  hilf <- DP_cands[which(DP_cands$v==v),]
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     nconst=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf)){
      k <- hilf[i,"k"]
      nconstk <- hilf$nconst[i]
      typek <- hilf[i,"type"]
      codek <- hilf[i,"code"]
    for (j in i:nrow(hilf)){
      l <- hilf[j,"k"]
      nconstl <- hilf$nconst[j]
      N <- hilf$N[i] + hilf$N[j] - min(v, nconstk + nconstl)
      typel <- hilf[j,"type"]
      codel <- hilf[j,"code"]
      kneu <- k*l
      codeneu <- paste0("productCA(",codek, ", ",
                        codel, ", c1=", nconstk, ")")
      nconstneu <- max(1, nconstk + nconstl - v)
        # in most cases this is accurate, in a few cases it can be slightly improved
        # but at high computational cost
        # length(attr(maxconstant(eval(parse(text=codeneu)), verbose=2), "constant_rows")$row_set_list)
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N,
                                     type=paste0(typek,"x",typel),
                                       nconst=nconstneu,
                                     code=codeneu,
                                     stage=2L))
    }
  }
  DPcat <- rbind(DPcat, addv)
}

rownames(DPcat) <- NULL
dim(DPcat)  ## 606 x 7

###### remove duplicated stage 2 cases #####
##  by keeping the one(s) with largest k and/or k1
## and the same or smaller N
uniques <- unique(DPcat[,c("v","N")])
dim(uniques)  ## 645x2
for (i in 1:nrow(uniques)){
  pick <- which(  DPcat$v==uniques$v[i] &
                    DPcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DPcat$k[pick])
  keepi1 <- which.max(DPcat$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DPcat <- DPcat[-remove,]
}
dim(DPcat) ## 455 rows
rownames(DPcat) <- NULL

######## remove dominated stage 2 rows, if any ####
remove <- numeric(0)
for (i in 1:(nrow(DPcat)-1)){
  ki <- DPcat$k[i]
  ci <- DPcat$nconst[i]

  vi <- DPcat$v[i]
  Ni <- DPcat$N[i]
  betteri <- which((DPcat[1:(i-1),]$k > ki |
                      DPcat[1:(i-1),]$nconst > ci |
                      DPcat[1:(i-1),]$N < Ni ) &
                     DPcat[1:(i-1),]$k >= ki &
                     DPcat[1:(i-1),]$nconst >= ci &
                     DPcat[1:(i-1),]$v == vi &
                     DPcat[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
DPcat <- DPcat[-remove,]
dim(DPcat) ## 452x7
rownames(DPcat) <- NULL
table(DPcat$stage)

# ###### check coverage for stage 2 ###############
# for (i in 1:nrow(DPcat)){
#   ## dimensions checked completely,
#   ## coverage only up to approx 335
#   ## should be OK for the others as well
#   ## all principles were now involved
#   print(i)
#   plan <- eval(parse(text=DPcat$code[i]))
#   if (!nrow(plan)==DPcat$N[i] && ncol(plan)==DPcat$k[i]){
#     print("wrong dimensions")
#     break
#   }
#   # if (!all(coverage(plan,2)==1)){
#   #   print("coverage violated")
#   #   break
#   # }
# }

##### add SCA_Bose entries to DP_cands ###############
## recursiveBose takes care of pairs of them, therefore only now
for (v in unique(DP_cands$v)){
  if (!v %in% primedat$q) next
  DP_cands <- rbind(DP_cands,
                     data.frame(v=v, k=v, N=v^2,
                                type="Bose",
                                nconst=v,
                                code=paste0("SCA_Bose(", v, ")[,1:", v, "]"),
                                stage=1L))
}
dim(DP_cands)  ## 140x7

###### remove duplicated stage 1 rows ######
##  by keeping the one(s) with largest k and/or nconst
uniques <- unique(DP_cands[,c("v","N")])
dim(uniques)  ## 136
for (i in 1:nrow(uniques)){
  pick <- which(  DP_cands$v==uniques$v[i] &
                    DP_cands$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DP_cands$k[pick])
  keepi1 <- which.max(DP_cands$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DP_cands <- DP_cands[-remove,]
}
dim(DP_cands) ## 141 rows from stage1, part 2

######## remove dominated stage 1 rows, if any ####
remove <- numeric(0)
for (i in 1:(nrow(DP_cands)-1)){
  ki <- DP_cands$k[i]
  ci <- DP_cands$nconst[i]

  vi <- DP_cands$v[i]
  Ni <- DP_cands$N[i]
  betteri <- which((DP_cands[(i+1):nrow(DP_cands),]$k > ki |
                      DP_cands[(i+1):nrow(DP_cands),]$nconst > ci |
                      DPcat[1:(i-1),]$N < Ni ) &
                     DP_cands[(i+1):nrow(DP_cands),]$k >= ki &
                     DP_cands[(i+1):nrow(DP_cands),]$nconst >= ci &
                     DP_cands[(i+1):nrow(DP_cands),]$v == vi &
                     DP_cands[(i+1):nrow(DP_cands),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
## nothing redundant
if (length(remove)>0)
  DP_cands <- DP_cands[-remove,]
dim(DP_cands) ## 139x7
rownames(DP_cands) <- NULL

# ############################ stage 3 ###################################
for (v in unique(DPcat$v)){
  ## loop over v
  hilf2 <- DPcat[which(DPcat$v==v & DPcat$N<=20000 & DPcat$stage==2),]  ## stage 2
  hilf1 <- DP_cands[which(DP_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0), type=character(0),
                     nconst=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    nconstk <- hilf1$nconst[i]
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      nconstl <- hilf2$nconst[j]
      N <- hilf1$N[i] + hilf2$N[j] - min(v, nconstk + nconstl)
      kneu <- k*l
      codeneu <- paste0("productCA(",codek, ", ",
                        codel, ", c1=", nconstk, ")")
      nconstneu <- max(1, nconstk + nconstl - v)
      # in most cases this is accurate, in a few cases it can be slightly improved
      # but at high computational cost
      # length(attr(maxconstant(eval(parse(text=codeneu)), verbose=2), "constant_rows")$row_set_list)
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N,
                                     type=paste0(typek,"x",typel),
                                     nconst=nconstneu,
                                     code=codeneu,
                                     stage=3L))
        }
  }
  DPcat <- rbind(DPcat, addv)
}
dim(DPcat) ## 4454 rows

###### redo the duplication removal #################
uniques <- unique(DPcat[,c("v","N")])
dim(uniques)  #660x2

for (i in 1:nrow(uniques)){
  pick <- which(  DPcat$v==uniques$v[i] &
                    DPcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DPcat$k[pick])
  keepi1 <- which.max(DPcat$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DPcat <- DPcat[-remove,]
}
dim(DPcat) ## 1168 rows, 405 from stage 2
table(DPcat$stage)
rownames(DPcat) <- NULL

###### check for dominated rows ###############################################
## remove dominated rows, if any
remove <- numeric(0)
for (i in 1:(nrow(DPcat)-1)){
  ki <- DPcat$k[i]
  ci <- DPcat$nconst[i]

  vi <- DPcat$v[i]
  Ni <- DPcat$N[i]
  betteri <- which((DPcat[1:(i-1),]$k > ki |
                      DPcat[1:(i-1),]$nconst > ci |
                      DPcat[1:(i-1),]$N < Ni ) &
                     DPcat[1:(i-1),]$k >= ki &
                     DPcat[1:(i-1),]$nconst >= ci &
                     DPcat[1:(i-1),]$v == vi &
                     DPcat[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}

DPcat <- DPcat[-remove,]
dim(DPcat) # 1141x7
## limit to k<=20000
## DPcat <- DPcat[DPcat$k<=20000,]
dim(DPcat) # 1045x7
rownames(DPcat) <- NULL

# ###################### stage 4 #############################
## add fourth factors among arrays with the same v
for (v in unique(DPcat$v)){
  ## loop over v
  hilf2 <- DPcat[which(DPcat$v==v & DPcat$N<=20000 & DPcat$stage==3),]  ## stage 3
  hilf1 <- DP_cands[which(DP_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0), type=character(0),
                     nconst=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    nconstk <- hilf1$nconst[i]
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      nconstl <- hilf2$nconst[j]
      N <- hilf1$N[i] + hilf2$N[j] - min(v, nconstk + nconstl)
      kneu <- k*l
      codeneu <- paste0("productCA(",codek, ", ",
                        codel, ", c1=", nconstk, ")")
      nconstneu <- max(1, nconstk + nconstl - v)
      # in most cases this is accurate, in a few cases it can be slightly improved
      # but at high computational cost
      # length(attr(maxconstant(eval(parse(text=codeneu)), verbose=2), "constant_rows")$row_set_list)
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N,
                                     type=paste0(typek,"x",typel),
                                     nconst=nconstneu,
                                     code=codeneu,
                                     stage=4L))
    }
  }
  DPcat <- rbind(DPcat, addv)
}
dim(DPcat) ## 6969 rows
table(DPcat$stage)
## DPcat <- DPcat[DPcat$k<=20000,]
dim(DPcat)  ## 3593 rows

###### redo the duplication removal within stage 4 ####################
uniques <- unique(DPcat[,c("v","N")])
dim(uniques)  # 740x2

for (i in 1:nrow(uniques)){
  pick <- which(  DPcat$v==uniques$v[i] &
                    DPcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DPcat$k[pick])
  keepi1 <- which.max(DPcat$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DPcat <- DPcat[-remove,]
}
dim(DPcat) ## 1342 rows
table(DPcat$stage)
rownames(DPcat) <- NULL

###### check for dominated rows ###############################################
## remove dominated rows, if any
remove <- numeric(0)
for (i in 1:(nrow(DPcat)-1)){
  ki <- DPcat$k[i]
  ci <- DPcat$nconst[i]

  vi <- DPcat$v[i]
  Ni <- DPcat$N[i]
  betteri <- which((DPcat[1:(i-1),]$k > ki |
                      DPcat[1:(i-1),]$nconst > ci |
                      DPcat[1:(i-1),]$N < Ni ) &
                     DPcat[1:(i-1),]$k >= ki &
                     DPcat[1:(i-1),]$nconst >= ci &
                     DPcat[1:(i-1),]$v == vi &
                     DPcat[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
DPcat <- DPcat[-remove,]
dim(DPcat) # 1307x7

lapply(unique(DPcat$v), function(obj) fivenum(DPcat$k[DPcat$v==obj]))

# ###################### stage 5 #############################
## add fifth factors among arrays with the same v
for (v in unique(DPcat$v)){
  ## loop over v
  hilf2 <- DPcat[which(DPcat$v==v & DPcat$N<=20000 & DPcat$stage==4),]
  hilf1 <- DP_cands[which(DP_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     nconst=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    nconstk <- hilf1$nconst[i]
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      nconstl <- hilf2$nconst[j]
      N <- hilf1$N[i] + hilf2$N[j] - min(v, nconstk + nconstl)
      kneu <- k*l
      codeneu <- paste0("productCA(",codek, ", ",
                        codel, ", c1=", nconstk, ")")
      nconstneu <- max(1, nconstk + nconstl - v)
      # in most cases this is accurate, in a few cases it can be slightly improved
      # but at high computational cost
      # length(attr(maxconstant(eval(parse(text=codeneu)), verbose=2), "constant_rows")$row_set_list)
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N,
                                     type=paste0(typek,"x",typel),
                                     nconst=nconstneu,
                                     code=codeneu,
                                     stage=5L))
    }
  }
  DPcat <- rbind(DPcat, addv)
}
dim(DPcat)  ## 4837 rows
table(DPcat$stage)
## DPcat <- DPcat[DPcat$k<=20000,]
dim(DPcat)  ## 2384 rows

###### redo the duplication removal within stage 5 ####################
uniques <- unique(DPcat[,c("v","N")])
dim(uniques)  # 739x2

for (i in 1:nrow(uniques)){
  pick <- which(  DPcat$v==uniques$v[i] &
                    DPcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DPcat$k[pick])
  keepi1 <- which.max(DPcat$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DPcat <- DPcat[-remove,]
}
dim(DPcat) ## 1397 rows
table(DPcat$stage)
rownames(DPcat) <- NULL

###### check for dominated rows ###############################################
## remove dominated rows, if any
remove <- numeric(0)
for (i in 1:(nrow(DPcat)-1)){
  ki <- DPcat$k[i]
  ci <- DPcat$nconst[i]

  vi <- DPcat$v[i]
  Ni <- DPcat$N[i]
  betteri <- which((DPcat[1:(i-1),]$k > ki |
                      DPcat[1:(i-1),]$nconst > ci |
                      DPcat[1:(i-1),]$N < Ni ) &
                     DPcat[1:(i-1),]$k >= ki &
                     DPcat[1:(i-1),]$nconst >= ci &
                     DPcat[1:(i-1),]$v == vi &
                     DPcat[1:(i-1),]$N <= Ni) + i
  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
DPcat <- DPcat[-remove,]
dim(DPcat) # 1347x7
rownames(DPcat) <- NULL

# ###################### stage 6 #############################
## add sixth factors among arrays with the same v
for (v in unique(DPcat$v)){
  ## loop over v
  hilf2 <- DPcat[which(DPcat$v==v & DPcat$N<=20000 & DPcat$stage==5),]
  hilf1 <- DP_cands[which(DP_cands$v==v),]
  if (nrow(hilf2)==0) next
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     type=character(0),
                     nconst=numeric(0), code=character(0),
                     stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    k <- hilf1[i,"k"]
    nconstk <- hilf1$nconst[i]
    typek <- hilf1$type[i]
    codek <- hilf1$code[i]
    for (j in 1:nrow(hilf2)){
      l <- hilf2[j,"k"]
      typel <- hilf2$type[j]
      codel <- hilf2$code[j]
      nconstl <- hilf2$nconst[j]
      N <- hilf1$N[i] + hilf2$N[j] - min(v, nconstk + nconstl)
      kneu <- k*l
      codeneu <- paste0("productCA(",codek, ", ",
                        codel, ", c1=", nconstk, ")")
      nconstneu <- max(1, nconstk + nconstl - v)
      # in most cases this is accurate, in a few cases it can be slightly improved
      # but at high computational cost
      # length(attr(maxconstant(eval(parse(text=codeneu)), verbose=2), "constant_rows")$row_set_list)
      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N,
                                     type=paste0(typek,"x",typel),
                                     nconst=nconstneu,
                                     code=codeneu,
                                     stage=6L))
    }
  }
  DPcat <- rbind(DPcat, addv)
}
dim(DPcat)  ## 2576 rows
table(DPcat$stage)
# DPcat <- DPcat[DPcat$k<=20000,]
dim(DPcat)  ## 1440 rows

###### redo the duplication removal within stage 6 ####################
uniques <- unique(DPcat[,c("v","N")])
dim(uniques)  # 761x2

for (i in 1:nrow(uniques)){
  pick <- which(  DPcat$v==uniques$v[i] &
                    DPcat$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(DPcat$k[pick])
  keepi1 <- which.max(DPcat$nconst[pick])
  ## remove all others
  remove <- pick[-c(keepi, keepi1)]
  ## remove the duplicated candidates
  if (length(remove) > 0)
    DPcat <- DPcat[-remove,]
}
dim(DPcat) ## 1353 rows, 28 from stage 6
table(DPcat$stage)
rownames(DPcat) <- NULL

###### final check for dominated rows #####################
###### nconst does not matter any more
## remove dominated rows, if any
remove <- numeric(0)
for (i in 1:nrow(DPcat)){
  ki <- DPcat$k[i]
  vi <- DPcat$v[i]
  Ni <- DPcat$N[i]
  betteri <- which((DPcat$k > ki |
                      DPcat$N < Ni) &
                     DPcat$v == vi &
                     DPcat$N <= Ni &
                      DPcat$k >=ki)

  if (length(betteri)>0){
    print(paste0("i:",i))
    print(betteri)
    remove <- c(remove,i)
  }
}
if (length(remove)>0) DPcat <- DPcat[-remove,]
dim(DPcat) # 528x7

DPcat <- DPcat[ord(DPcat[,c("v","k")]),]
rownames(DPcat) <- NULL

sapply(unique(DPcat$v), function(obj) fivenum(DPcat$k[DPcat$v==obj]))
sapply(unique(DPcat$v), function(obj) fivenum(DPcat$N[DPcat$v==obj]))

## limit so that k=20000 remains covered
remove <- numeric(0)
for (v in unique(DPcat$v)){
  hilf <- which(DPcat$v==v & DPcat$k > 20000)
  remove <- c(remove, setdiff(hilf, hilf[which.min(DPcat[hilf,]$N)]))
}
DPcat <- DPcat[-remove,]
dim(DPcat)
rownames(DPcat) <- NULL

table(isbest <- bestNvec(DPcat$k, DPcat$v)==DPcat$N, DPcat$v)
table(isbest, DPcat$stage)
table(isopt <- eCANvec(DPcat$k, DPcat$v)==DPcat$N, DPcat$v)
table(isopt, DPcat$stage)

######## check dimensions for everything ################
## stage
table(DPcat$stage)
for (i in 1:nrow(DPcat)){
  print(i)
  plan <- dpCA(DPcat$k[i], DPcat$v[i])
  if (!nrow(plan)==DPcat$N[i] && ncol(plan)==DPcat$k[i]){
    if (nrow(plan)<DPcat$N[i]) message(DPcat$N[i]-nrow(plan), " fewer rows") else{
    print("wrong dimensions")
    break
    }
  }
}

DPcat[c(49, 59, 62, 86, 97, 109),]$N <- DPcat[c(49, 59, 62, 86, 97, 109),]$N - 1
## in rows 49, 59, 62, 86, 97, 109 N was one less than claimed
DPcat[c(114, 115, 350),]$N <- DPcat[c(114, 115, 350),]$N - 2
## in rows 114, 115, 350,  N was two less than claimed


boxplot(DPcat$N/eCANvec(DPcat$k, DPcat$v)~DPcat$v, las=1, horizontal=TRUE)
boxplot(DPcat$N/bestNvec(DPcat$k, DPcat$v)~DPcat$v, las=1, horizontal=TRUE)

# save(DP_cands, file="D:/rtests/CAs/inst/extdata/DP_cands.rda", compress="xz")
# save(DPcat, file="D:/rtests/CAs/data/DPcat.rda", compress="xz")
