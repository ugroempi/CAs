## !!! NOTE: FOR REPRODUCING CREATION OF THE FILES,
## !!! the production function CS_CMMSSY must be modified
## !!! to forget restricting the files to k columns (esp. in
## !!!                the return command)

## typed in from examples of Colbourn et al. 2006 (CMMSSY)
##   start with Inf
## to be used with the Z_v-1
## first level: v
## second level: k-1 (length of starter)
CMMSSYStarters <- list(
  '5' = list( '5'=c(Inf,0,1,3,2),
              '6'=c(Inf,0,0,1,3,2),
              '7'=c(Inf,0,0,0,1,3,2),
              '8'=c(Inf,0,0,0,0,1,3,2),
              '9'=c(Inf,0,0,0,0,0,1,3,2),
             '10'=c(Inf,0,0,0,0,0,0,1,3,2)),
  '6' = list( '8'=c(Inf,0,1,3,0,2,1,4),
              '9'=c(Inf,0,0,1,0,0,3,2,4),
             '10'=c(Inf,0,0,0,1,1,4,3,0,2),
             '11'=c(Inf,0,0,0,0,1,0,4,2,3,0),
             '12'=c(Inf,0,0,0,0,0,1,0,4,2,3,0)),
  '7' = list( '7'=c(Inf,0,2,1,4,5,3),
              '9'=c(Inf,0,0,2,1,4,5,3,3),
             '10'=c(Inf,0,0,0,1,0,3,5,4,2),
             '11'=c(Inf,0,0,0,0,1,0,3,5,4,2),
             '12'=c(Inf,0,0,0,0,0,1,0,3,5,4,2),
             '13'=c(Inf,0,0,0,0,0,0,1,0,3,5,4,2)),
  '8' = list('11'=c(Inf,0,1,0,0,0,3,5,4,2,6),
             '12'=c(Inf,0,0,0,1,0,6,4,6,3,2,5),
             '13'=c(Inf,0,0,0,0,1,2,6,5,3,6,2,4),
             '14'=c(Inf,0,0,0,0,0,1,0,4,6,4,3,2,5)),
  '9' = list('12'=c(Inf,0,1,6,4,5,0,7,6,2,1,3),
             '13'=c(Inf,0,0,0,0,2,6,1,7,6,3,4,5),
             '14'=c(Inf,0,0,0,0,0,2,6,1,7,6,3,4,5),
             '15'=c(Inf,0,0,0,0,0,0,2,6,1,7,6,3,4,5))
)

replacefun <- function(tobereplaced, replacement, tobetreated,
                       Nj, kj, k1j, Nkeepi, kkeepi, k1keepi
                       ){
  ## function to replace one of the two ingredients in a productPCA call
  ## tobereplaced is the code entry that is to be changed
  ## replacement is the code entry it should be changed to
  ## tobetreated is a data frame whose code entry must be treated
  ##    for the rows in which the tobereplaced is found
  ## ... j entries refer to the to be replaced prior status (numbers)
  ## ... keepi entries refer to the replacement status (numbers)

  ## !!! the function assumes that this is done at each step, so that
  ##     hopefully the problems always refer to the highest level application
  ##     of productPCA
  gefunden <- grep(tobereplaced, tobetreated$code, fixed=TRUE)
  gefunden2 <- grep(tobereplaced, sub(tobereplaced, "xxx", tobetreated$code, fixed=TRUE), fixed=TRUE)
  if (length(gefunden)==0) return(tobetreated)
  stopifnot(kkeepi > kj)
  stopifnot(Nj<=Nkeepi)
  tobetreated$code[gefunden] <- gsub(pattern=tobereplaced, replacement = replacement,
                                      x = tobetreated$code[gefunden], fixed=TRUE)
  kalt <- tobetreated$k[gefunden]
  Nalt <- tobetreated$N[gefunden]
  k1alt <- tobetreated$k1[gefunden]

  l1 <- k1alt / k1j ## calculation
  ## k = k1*l1 + k2*l1 + k1*l2 = k*l1 + k1*l2
  l2 <- (kalt - kj*l1)/k1j
  kneu <- kkeepi*l1 + k1keepi*l2
  k1neu <- k1keepi*l1
  Nneu <- Nalt - Nj + Nkeepi
  tobetreated$k[gefunden] <- kneu
  tobetreated$k1[gefunden] <- k1neu
  tobetreated$N[gefunden] <- Nneu

  ## handle case where tobereplaced is contained twice
  if (length(gefunden2)>0){
    ## the entry occurs twice (possible only at stage 2)
    kalt <- tobetreated$k[gefunden2]
    Nalt <- tobetreated$N[gefunden2]
    valt <- tobetreated$v[gefunden2]
    k1alt <- tobetreated$k1[gefunden2]

    l1 <- k1alt / k1j ## calculation
    ## k = k1*l1 + k2*l1 + k1*l2 = k*l1 + k1*l2
    l2 <- (kalt - kj*l1)/k1j
    kneu <- kkeepi*l1 + k1keepi*l2
    k1neu <- k1keepi*l1
    Nneu <- Nalt - Nj + Nkeepi
    tobetreated$k[gefunden2] <- kneu
    tobetreated$k1[gefunden2] <- k1neu
    tobetreated$N[gefunden2] <- Nneu
  }
  tobetreated
}

##### initialize CMMSSYCombis with starter constructions (stage 1) #########
CMMSSYCombis <- cbind(
  v=(v <- rep(as.numeric(names(CMMSSYStarters)), lengths(CMMSSYStarters))),
  k=(k <- unlist(lapply(CMMSSYStarters, function(obj) as.numeric(names(obj))))+1),
  N=(v-1)*(k-1) + v,
  k1=(k1 <- k-1)
)
CMMSSYCombis <- as.data.frame(CMMSSYCombis)
## add code column for the starter rows
CMMSSYCombis <- cbind(CMMSSYCombis, code=paste0("CS_CMMSSY(", k, ",", v, ", starter=CMMSSYStarters[[as.character(", v,")]][[as.character(",k,"-1)]])"),
                      stage=1L)

######### stage 2 ####################
## add pairwise combinations among arrays with the same v
for (v in unique(CMMSSYCombis$v)){
hilf <- CMMSSYCombis[which(CMMSSYCombis$v==v),]
addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                   k1=numeric(0), code=character(0), stage=numeric(0))
for (i in 1:nrow(hilf)){
  for (j in i:nrow(hilf)){
    k <- hilf[i,"k"]
    l <- hilf[j,"k"]
    N <- hilf$N[i] + hilf$N[j] - v
    k1 <- hilf$k1[i]
    l1 <- hilf$k1[j]
    k2 <- k - k1
    l2 <- l - l1
    kneu <- k1*l1 + k1*l2 + k2*l1
    k1neu <- k1*l1
    addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                          code=paste0("productPCA(CS_CMMSSY(", k, "," , v,
                                      "), CS_CMMSSY(",l, ",", v, "))"),
                          stage=2L))
}
}
CMMSSYCombis <- rbind(CMMSSYCombis, addv)
}

rownames(CMMSSYCombis) <- NULL
dim(CMMSSYCombis)  ## 102 x 6

###### remove duplicates and dominated cases by keeping the one with largest k #####
## and the same or smaller N
uniques <- unique(CMMSSYCombis[,c("v","N")])
dim(uniques)
for (i in 1:nrow(uniques)){
  pick <- which(  CMMSSYCombis$v==uniques$v[i] &
                  CMMSSYCombis$N==uniques$N[i])
  if (length(pick)==1) next

  ## caution: keepi refers to pick as its reference set
  keepi <- which.max(CMMSSYCombis$k[pick])
  ## remove all others
  remove <- pick[-keepi]
  ## entries that can occur in the code column for other entries
  removestage1 <- intersect(remove, which(CMMSSYCombis$stage==1))

  if (length(removestage1) > 0){
    ## replace the run sizes and code for entries in stage 2 that used the removed stage 1 entries
    ## in order to avoid later disconnects
    for (j in 1:length(removestage1)){
      kj <- CMMSSYCombis$k[removestage1[j]]
      vj <- CMMSSYCombis$v[removestage1[j]]
      Nj <- CMMSSYCombis$N[removestage1[j]]
      k1j <- CMMSSYCombis$k1[removestage1[j]]

      tobereplaced <- paste0("CS_CMMSSY(", kj, "," , vj, ")")

      kkeepi <- CMMSSYCombis$k[pick[keepi]]
      Nkeepi <- CMMSSYCombis$N[pick[keepi]]
      vkeepi <- CMMSSYCombis$v[pick[keepi]]
      k1keepi <- CMMSSYCombis$k1[pick[keepi]]
      stopifnot(vj==vkeepi)
      replacement <- paste0("CS_CMMSSY(", kkeepi, "," , vkeepi, ")")
      ## handle the entries of CMMSSYCombis that contain the current removal
      ## candidate
      CMMSSYCombis <- replacefun(tobereplaced, replacement, tobetreated=CMMSSYCombis,
                                 Nj, kj, k1j, Nkeepi, kkeepi, k1keepi)
    }
  }
  ## remove the duplicated candidates
  if (length(remove) > 0)
      CMMSSYCombis <- CMMSSYCombis[-remove,]
}
dim(CMMSSYCombis) ## 70 rows, 24 from stage1

## checked, coverage and all dimensions are correct up to here

## remove dominated rows, if any (none at stage 2 after removal of duplicated stage 1)
for (i in 1:(nrow(CMMSSYCombis)-1)){
  ki <- CMMSSYCombis$k[i]
  vi <- CMMSSYCombis$v[i]
  Ni <- CMMSSYCombis$N[i]
  betteri <- which(CMMSSYCombis[(i+1):nrow(CMMSSYCombis),]$k > ki &
                     CMMSSYCombis[(i+1):nrow(CMMSSYCombis),]$v == vi &
                     CMMSSYCombis[(i+1):nrow(CMMSSYCombis),]$N <= Ni)
  if (length(betteri)>0) print(i)
}

###### check coverage for up to stage 2 ###############
# MSCs <- CMMSSYCombis
# the checks below ran without throwing an error
# for up to stage 2

# vs <- MSCs$v; ks <- MSCs$k; Ns <- MSCs$N
#
# for (i in 1:nrow(MSCs)){
#   print(i)
#   hilf <- CS_CMMSSY(ks[i], vs[i])
#   if (!nrow(hilf)==Ns[i] && ncol(hilf)==ks[i]) stop("Wrong dimensions for row ", i)
#   if (!all(coverage(hilf, 2, parallel=ifelse(ks[i]<=100, 1, 4))==1)) stop("Wrong code entry for row ", i)
# }

# ############################ stage 3 ###################################
for (v in unique(CMMSSYCombis$v)){
  ## loop over v (5 to 9)
  hilf <- CMMSSYCombis[which(CMMSSYCombis$v==v),]
  hilf1 <- hilf[which(hilf$stage==1),]  ## arrays with a starter, 24 rows
  hilf2 <- hilf[which(hilf$stage==2),]  ## arrays from stage 2
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     k1=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    for (j in 1:nrow(hilf2)){
      k <- hilf1[i,"k"]
      l <- hilf2[j,"k"]
      N <- hilf1$N[i] + hilf2$N[j] - v
      k1 <- hilf1$k1[i]
      l1 <- hilf2$k1[j]
      k2 <- k - k1
      l2 <- l - l1
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1

      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                                     code=paste0("productPCA(",
                                     ifelse(length(grep("starter",hilf1$code[i],fixed=TRUE))>0,
                                            paste0("CS_CMMSSY(", k, "," , v,
                                                 ")"), hilf1$code[i]),
                                     ", ",
                                     ifelse(length(grep("starter",hilf2$code[j],fixed=TRUE))>0,
                                            paste0("CS_CMMSSY(",l, ",", v, ")"),
                                            hilf2$code[j]),
                                     ")"),
                                     stage=3L))
    }
    }
  CMMSSYCombis <- rbind(CMMSSYCombis, addv)
}
dim(CMMSSYCombis) ## 298 rows

###### redo the duplication removal #################
uniques <- unique(CMMSSYCombis[,c("v","N")])
dim(uniques)
for (i in 1:nrow(uniques)){
  pick <- which(  CMMSSYCombis$v==uniques$v[i] &
                  CMMSSYCombis$N==uniques$N[i])
  if (length(pick)==1) next
  keepi <- which.max(CMMSSYCombis$k[pick])
  keepcode <- CMMSSYCombis$code[pick[keepi]]
  remove <- pick[-keepi]
  ## cases which could occur in construction (assuming stage 1 is not dominated by stage 3)
  removestage1 <- intersect(remove, which(CMMSSYCombis$stage==2))

  if (length(removestage1) > 0){
    ## replace the run sizes and code for entries in stage 1 that are used in later stage
    ## in order to avoid later disconnects
    for (j in 1:length(removestage1)){
      kj <- CMMSSYCombis$k[removestage1[j]]
      vj <- CMMSSYCombis$v[removestage1[j]]
      Nj <- CMMSSYCombis$N[removestage1[j]]
      k1j <- CMMSSYCombis$k1[removestage1[j]]

      tobereplaced <- CMMSSYCombis$code[removestage1[j]]

      kkeepi <- CMMSSYCombis$k[pick[keepi]]
      Nkeepi <- CMMSSYCombis$N[pick[keepi]]
      vkeepi <- CMMSSYCombis$v[pick[keepi]]
      k1keepi <- CMMSSYCombis$k1[pick[keepi]]
      stopifnot(vj==vkeepi)
      replacement <- CMMSSYCombis$code[pick[keepi]]
  print(rbind(tobereplaced=tobereplaced,
              replacement=replacement))
      ## handle the entries of CMMSSYCombis that contain the current removal
      ## candidate
      CMMSSYCombis <- replacefun(tobereplaced, replacement, tobetreated=CMMSSYCombis,
                                 Nj, kj, k1j, Nkeepi, kkeepi, k1keepi)
    }
  }
  ## remove the duplicated candidates
 # if (length(removestage1)>0) print(CMMSSYCombis[remove,])
  if (length(remove) > 0)
    CMMSSYCombis <- CMMSSYCombis[-remove,]
}
dim(CMMSSYCombis) ## 124 rows, 24 from stage1
rownames(CMMSSYCombis) <- NULL

## dominated cases
MSCs <- CMMSSYCombis
vs <- MSCs[,1]; ks <- MSCs[,2]; Ns <- MSCs[,3]; code <- MSCs$code
for (i in 1:nrow(MSCs)){
  hilf <- CS_CMMSSY(k=ks[i],v=vs[i])
  fromi <- which(CMMSSYCombis$k==ncol(hilf) & CMMSSYCombis$v==vs[i] & CMMSSYCombis$N==nrow(hilf))
  if (!(nrow(hilf)==Ns[i] && ncol(hilf)==ks[i])){
    print(dim(hilf))
    print(rbind(MSCs[i,],
                CMMSSYCombis[fromi,]))
  }
}

###### check for dominated rows ###############################################
## remove dominated rows, if any
## at stage 3: 46 is worse than 87
remove <- numeric(0)
removals <- vector(mode="list")
for (i in 1:(nrow(CMMSSYCombis)-1)){
  ki <- CMMSSYCombis$k[i]
  vi <- CMMSSYCombis$v[i]
  Ni <- CMMSSYCombis$N[i]
  betteri <- i + which(CMMSSYCombis[(i+1):nrow(CMMSSYCombis),]$k > ki &
                         CMMSSYCombis[(i+1):nrow(CMMSSYCombis),]$v == vi &
                         CMMSSYCombis[(i+1):nrow(CMMSSYCombis),]$N <= Ni)
  if (length(betteri)>0){
    print(i)
    print(betteri)
    remove <- c(remove,i)
    removals[[length(removals)+1]] <- list(better=CMMSSYCombis$code[betteri],
                                           removed=CMMSSYCombis$code[i])
  }
}
if (length(remove)>0){
  saferemove <- numeric(0)
  for (i in remove){
    if (all(grep(CMMSSYCombis$code[i], CMMSSYCombis$code, fixed=TRUE)==i))
      saferemove <- c(saferemove, i)
  }
  CMMSSYCombis <- CMMSSYCombis[-saferemove,]
}

dim(CMMSSYCombis)  ## 123 rows
rownames(CMMSSYCombis) <- NULL

## dominated in reality
MSCs <- CMMSSYCombis
vs <- MSCs[,1]; ks <- MSCs[,2]; Ns <- MSCs[,3]; code <- MSCs$code
for (i in 1:nrow(MSCs)){
  hilf <- CS_CMMSSY(k=ks[i],v=vs[i])
  fromi <- which(CMMSSYCombis$k==ncol(hilf) & CMMSSYCombis$v==vs[i] & CMMSSYCombis$N==nrow(hilf))
  if (!(nrow(hilf)==Ns[i] && ncol(hilf)==ks[i])){
    print(dim(hilf))
    print(rbind(MSCs[i,],
                CMMSSYCombis[fromi,]))
  }
}
## row 97 is dominated by row 91, both stage 3 --> remove
CMMSSYCombis <- CMMSSYCombis[-97,]
dim(CMMSSYCombis) ## 122 rows

# ###################### stage 4 #############################
## add fourth factors among arrays with the same v
for (v in unique(CMMSSYCombis$v)){
  ## loop over v (5 to 9)
  hilf <- CMMSSYCombis[which(CMMSSYCombis$v==v),]
  hilf1 <- hilf[which(hilf$stage==1),]
  hilf2 <- hilf[which(hilf$stage==3),]
  addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
                     k1=numeric(0), code=character(0), stage=numeric(0))
  for (i in 1:nrow(hilf1)){
    for (j in 1:nrow(hilf2)){
      k <- hilf1[i,"k"]
      l <- hilf2[j,"k"]
      N <- hilf1$N[i] + hilf2$N[j] - v
      k1 <- hilf1$k1[i]
      l1 <- hilf2$k1[j]
      k2 <- k - k1
      l2 <- l - l1
      kneu <- k1*l1 + k1*l2 + k2*l1
      k1neu <- k1*l1

      addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
                            code=paste0("productPCA(",
                                 ifelse(length(grep("starter",hilf1$code[i],fixed=TRUE))>0,
                                     paste0("CS_CMMSSY(", k, "," , v,
                                      ")"), hilf1$code[i]),
                                                 ", ",
                                  ifelse(length(grep("starter",hilf2$code[j],fixed=TRUE))>0,
                                      paste0("CS_CMMSSY(",l, ",", v, ")"),
                                        hilf2$code[j]),
                                       ")"),
                                     stage=4L))
    }
  }
  CMMSSYCombis <- rbind(CMMSSYCombis, addv)
}
dim(CMMSSYCombis)  ## 429 rows

###### redo the duplication removal within stage 4 ####################
uniques <- unique(CMMSSYCombis[,c("v","N")])
dim(uniques)
for (i in 1:nrow(uniques)){
  pick <- which(  CMMSSYCombis$v==uniques$v[i] &
                  CMMSSYCombis$N==uniques$N[i])
  if (length(pick)==1) next
  keepi <- which.max(CMMSSYCombis$k[pick])
  remove <- pick[-keepi]

  ## cases which could occur in construction (assuming stage 1 and 2 are not dominated by stage 4
  removestage1 <- intersect(remove, which(CMMSSYCombis$stage==3))

  if (length(removestage1) > 0){
    ## replace the run sizes and code for entries in stage 1 that are used in later stage
    ## in order to avoid later disconnects
    for (j in 1:length(removestage1)){
      kj <- CMMSSYCombis$k[removestage1[j]]
      vj <- CMMSSYCombis$v[removestage1[j]]
      Nj <- CMMSSYCombis$N[removestage1[j]]
      k1j <- CMMSSYCombis$k1[removestage1[j]]

      tobereplaced <- CMMSSYCombis$code[removestage1[j]]

      kkeepi <- CMMSSYCombis$k[pick[keepi]]
      Nkeepi <- CMMSSYCombis$N[pick[keepi]]
      vkeepi <- CMMSSYCombis$v[pick[keepi]]
      k1keepi <- CMMSSYCombis$k1[pick[keepi]]
      stopifnot(vj==vkeepi)
      replacement <- CMMSSYCombis$code[pick[keepi]]
print(rbind(tobereplaced=tobereplaced, replacement=replacement))
      ## handle the entries of CMMSSYCombis that contain the current removal
      ## candidate
      CMMSSYCombis <- replacefun(tobereplaced, replacement, tobetreated=CMMSSYCombis,
                                 Nj, kj, k1j, Nkeepi, kkeepi, k1keepi)
    }
  }
  ## remove the duplicated candidates
 # if (length(removestage1)>0) print(CMMSSYCombis[remove,])
  if (length(remove) > 0)
    CMMSSYCombis <- CMMSSYCombis[-remove,]
}
dim(CMMSSYCombis) ## 182 rows, 24 from stage1
rownames(CMMSSYCombis) <- NULL

###### dominated cases ##############################
# MSCs <- CMMSSYCombis
# vs <- MSCs[,1]; ks <- MSCs[,2]; Ns <- MSCs[,3]; code <- MSCs$code
# for (i in 1:nrow(MSCs)){
#   hilf <- CS_CMMSSY(k=ks[i],v=vs[i])
#   fromi <- which(CMMSSYCombis$k==ncol(hilf) & CMMSSYCombis$v==vs[i] & CMMSSYCombis$N==nrow(hilf))
#   if (!(nrow(hilf)==Ns[i] && ncol(hilf)==ks[i])){
#     print(dim(hilf))
#     print(rbind(MSCs[i,],
#                 CMMSSYCombis[fromi,]))
#   }
# }

## row 146 is dominated by row 145 (both stage 4)
## row 150 is dominated by row 139 (both stage 4)

CMMSSYCombis <- CMMSSYCombis[-c(146,150),] ## recursiveBose analogues from more smaller are better than from fewer larger
dim(CMMSSYCombis)

# ############## stage 5 is omitted ############
# as it yields very very large arrays
# for (v in unique(CMMSSYCombis$v)){
#   hilf <- CMMSSYCombis[which(CMMSSYCombis$v==v),]
#   hilf1 <- hilf[which(hilf$stage==1),]
#   hilf2 <- hilf[which(hilf$stage==4),]
#   addv <- data.frame(v=numeric(0), k=numeric(0), N=numeric(0),
#                      k1=numeric(0), code=character(0), stage=numeric(0))
#   for (i in 1:nrow(hilf1)){
#     for (j in 1:nrow(hilf2)){
#       k <- hilf1[i,"k"]
#       l <- hilf2[j,"k"]
#       N <- hilf1$N[i] + hilf2$N[j] - v
#       k1 <- hilf1$k1[i]
#       l1 <- hilf2$k1[j]
#       k2 <- k - k1
#       l2 <- l - l1
#       kneu <- k1*l1 + k1*l2 + k2*l1
#       k1neu <- k1*l1
#
#       addv <- rbind(addv, data.frame(v=v, k=kneu, N=N, k1=k1neu,
#                                code=paste0("productPCA(",
#                                            ifelse(length(grep("starter",hilf1$code[i],fixed=TRUE))>0,
#                                                   paste0("CS_CMMSSY(", k, "," , v,
#                                                          ")"), hilf1$code[i]),
#                                            ", ",
#                                            ifelse(length(grep("starter",hilf2$code[j],fixed=TRUE))>0,
#                                                   paste0("CS_CMMSSY(",l, ",", v, ")"),
#                                                   hilf2$code[j]),
#                                            ")"),
#                                stage=5L))
#  }
# }
#   CMMSSYCombis <- rbind(CMMSSYCombis, addv)
# }
# dim(CMMSSYCombis)  ## 515 rows
# rownames(CMMSSYCombis) <- NULL
#
# ## redo the duplication removal within stage 5
# uniques <- unique(CMMSSYCombis[,c("v","N")])
# dim(uniques)
# for (i in 1:nrow(uniques)){
#   pick <- which(  CMMSSYCombis$v==uniques$v[i] &
#                   CMMSSYCombis$N==uniques$N[i])
#   if (length(pick)==1) next
#   keepi <- which.max(CMMSSYCombis$k[pick])
#   remove <- pick[-keepi]
#
#   ## cases which could occur in construction (assuming stage 1 to 3 are not dominated by stage 5
#   removestage1 <- intersect(remove, which(CMMSSYCombis$stage==4))
#
#   if (length(removestage1) > 0){
#     ## replace the run sizes and code for entries in stage 1 that are used in later stage
#     ## in order to avoid later disconnects
#     for (j in 1:length(removestage1)){
#       kj <- CMMSSYCombis$k[removestage1[j]]
#       vj <- CMMSSYCombis$v[removestage1[j]]
#       Nj <- CMMSSYCombis$N[removestage1[j]]
#       k1j <- CMMSSYCombis$k1[removestage1[j]]
#
#       tobereplaced <- CMMSSYCombis$code[removestage1[j]]
#
#       kkeepi <- CMMSSYCombis$k[pick[keepi]]
#       Nkeepi <- CMMSSYCombis$N[pick[keepi]]
#       vkeepi <- CMMSSYCombis$v[pick[keepi]]
#       k1keepi <- CMMSSYCombis$k1[pick[keepi]]
#       stopifnot(vj==vkeepi)
#       replacement <- CMMSSYCombis$code[pick[keepi]]
#       print(rbind(tobereplaced=tobereplaced, replacement=replacement))
#       ## handle the entries of CMMSSYCombis that contain the current removal
#       ## candidate
#       CMMSSYCombis <- replacefun(tobereplaced, replacement, tobetreated=CMMSSYCombis,
#                                  Nj, kj, k1j, Nkeepi, kkeepi, k1keepi)
#     }
#   }
#   ## remove the duplicated candidates
#   # if (length(removestage1)>0) print(CMMSSYCombis[remove,])
#   if (length(remove) > 0)
#     CMMSSYCombis <- CMMSSYCombis[-remove,]
# }
# dim(CMMSSYCombis) ## 243 rows, 24 from stage1
#
# rownames(CMMSSYCombis) <- NULL
# table(CMMSSYCombis$stage)


######## check dimensions for all stages ################
## stage
MSCs <- CMMSSYCombis#[-c(44,72,73),]
vs <- MSCs[,1]; ks <- MSCs[,2]; Ns <- MSCs[,3]
for (i in 1:nrow(MSCs)){
  hilf <- CS_CMMSSY(k=ks[i],v=vs[i])
  if (!(nrow(hilf)==Ns[i] && ncol(hilf)==ks[i])) stop("Wrong dimensions for row ", i)
}

## dominated cases
MSCs <- CMMSSYCombis
rownames(MSCs) <- NULL
vs <- MSCs[,1]; ks <- MSCs[,2]; Ns <- MSCs[,3]; code <- MSCs$code
for (i in 1:nrow(MSCs)){
  hilf <- CS_CMMSSY(k=ks[i],v=vs[i])
  fromi <- which(CMMSSYCombis$k==ncol(hilf) & CMMSSYCombis$v==vs[i] & CMMSSYCombis$N==nrow(hilf))
  if (!(nrow(hilf)==Ns[i] && ncol(hilf)==ks[i])){
    print(dim(hilf))
    print(rbind(MSCs[i,],
                CMMSSYCombis[fromi,]))
  }
}
