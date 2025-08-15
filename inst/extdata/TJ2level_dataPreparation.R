## code for reading the 2-level files available in the extdata folder
## and storing them in TJ2level_CAs

## Author: Ulrike Groemping

# original location of the files in Feb 2025
# https://www.tamps.cinvestav.mx/~oc/CA/v2/t3/N54k968v2%5E968t3.ca
# source of the files: Jose Torres-Jimenez homepage

pfad_extdata <- system.file("extdata", package="CAs")
fileloc <- paste(pfad_extdata, "TJ2level", sep="/")
toberead <- list.files(fileloc)
fns <- rownames(TJcat)[which(TJcat[,"v"]==2 & !TJcat$nameInTJ2level_CAs=="")]
length(setdiff(toberead, fns))  ## should be 0
length(setdiff(fns, toberead))  ## should be 0

TJ2level_CAs <- vector(mode="list", length=length(fns))
m <- 0
for (i in 1:length(fns)){
    m <- m+1
    print(fns[i])
    TJ2level_CAs[[m]] <- readCA(paste0(fileloc,"/", fns[i]), flexible.symbols = 2,
                                origin="Torres-Jimenez repository, February 6 2025")
}
names(TJ2level_CAs) <- TJcat[fns,]$nameInTJ2level_CAs

# save(TJ2level_CAs, file="D:/rtests/CAs/data/TJ2level_CAs.rda", compress="xz")

## check that it works
tjCA(3,21,2)  ## construction from the package
tjCA(6,18,2)  ## stored array

### code that was used for removing arrays for reduced package size
### to be re-run sporadically for identifying new opportunities

### reduce the storage requirements, as these files are responsible for
###    more than 80pct of the package size

## handle rows without stored arrays
settings <- TJcat[which(TJcat$nameInTJ2level_CAs=="" & TJcat$replaceable=="" & TJcat$v==2),]
for (i in 1:nrow(settings)){
  hilf <- names(bestN(settings$t[i], settings$k[i], 2))
  if (!hilf=="TJ") settings$replaceable[i] <- hilf
}
settings <- settings[which(!settings$replaceable==""),] ## August 15 2025: 33 rows
                                                        ## CS_CK, CYCLOTOMY, PALEY, powerCT
TJcat[rownames(settings),]$replaceable <- settings$replaceable

## handle rows with stored arrays
settings <- TJcat[which(!TJcat$nameInTJ2level_CAs==""),]
settings$replaceable <- ""
for (i in 1:nrow(settings)){
  hilf <- names(bestN(settings$t[i], settings$k[i], 2))
  if (!hilf=="TJ") settings$replaceable[i] <- hilf
}

table(settings$replaceable)
## do not remove arrays that are only on the Internet
settings$replaceable[which(settings$replaceable=="DWYER")] <- ""
table(settings$replaceable) ## seven additional arrays removable July 10 2025
                            ## three additional arrays removable August 15 2025

## the following statement was only needed at the first round
# TJcat$replaceable <- ""
## August 15 2025
## modified replaceable for CA(12, 3, 11, 2) to miscCA, as this has two constant rows
TJcat["ca.3.2^11.txt",]$replaceable <- "miscCA"

TJcat[rownames(settings),]$replaceable <- settings$replaceable
## the following statement was only needed at the first round
# TJcat$code <- ""
TJcat[rownames(settings),]$code <- ""
## the following statement was only needed at the first round
# TJcat$code[TJcat$v==3 & TJcat$t==2] <- paste0("CAEX(", TJcat$k[TJcat$v==3 & TJcat$t==2], ")")

## powerCT CK_doublingCA CKRS CS_CK PALEY WKS
TJcat$code[TJcat$replaceable=="powerCT"] <- paste0("powerCA(",
                                                   TJcat$t[TJcat$replaceable=="powerCT"],
                                                   ", ",
                                                   TJcat$k[TJcat$replaceable=="powerCT"],
                                                   ", ",
                                                   TJcat$v[TJcat$replaceable=="powerCT"],
                                                         ")")
TJcat$code[TJcat$replaceable=="CK_doublingCA"] <- paste0("CK_doublingCA(",
                                                         TJcat$k[TJcat$replaceable=="CK_doublingCA"],
                                                         ", ",
                                                         TJcat$v[TJcat$replaceable=="CK_doublingCA"],
                                                         ")")
TJcat$code[TJcat$replaceable=="miscCA"] <- paste0("miscCA(",
                                                TJcat$t[TJcat$replaceable=="miscCA"],
                                                ", ",
                                                TJcat$k[TJcat$replaceable=="miscCA"],
                                                ", ",
                                                TJcat$v[TJcat$replaceable=="miscCA"],
                                                ")")
TJcat$code[TJcat$replaceable=="CYCLOTOMY"] <- paste0("cyclotomyCA(",
                                                TJcat$t[TJcat$replaceable=="CYCLOTOMY"],
                                                ", ",
                                                TJcat$k[TJcat$replaceable=="CYCLOTOMY"],
                                                ", ",
                                                TJcat$v[TJcat$replaceable=="CYCLOTOMY"],
                                                ")")

TJcat$code[TJcat$replaceable=="CKRS"] <- paste0("ckrsCA(",
                                                TJcat$t[TJcat$replaceable=="CKRS"],
                                                ", ",
                                                TJcat$k[TJcat$replaceable=="CKRS"],
                                                ", ",
                                                TJcat$v[TJcat$replaceable=="CKRS"],
                                                         ")")

TJcat$code[TJcat$replaceable=="CS_CK"] <- paste0("CS_CK(",
                                                TJcat$k[TJcat$replaceable=="CS_CK"],
                                                ", ",
                                                TJcat$t[TJcat$replaceable=="CS_CK"],
                                                ", ",
                                                TJcat$v[TJcat$replaceable=="CS_CK"],
                                                ")")

TJcat$code[TJcat$replaceable=="PALEY"] <- paste0("paleyCA(",
                                                 TJcat$t[TJcat$replaceable=="PALEY"],
                                                 ", ",
                                                 TJcat$k[TJcat$replaceable=="PALEY"],
                                                 ")")
TJcat$code[TJcat$replaceable=="WKS"] <- paste0('WKS_CAs[["', TJcat$k[TJcat$replaceable=="WKS"], '"]]')

TJcat$code[TJcat$replaceable=="" & TJcat$code=="" & !TJcat$nameInTJ2level_CAs==""] <-
  paste0("TJ2level_CAs[['", TJcat$nameInTJ2level_CAs[TJcat$replaceable=="" &
                           TJcat$code=="" & !TJcat$nameInTJ2level_CAs==""], "']]")
## files to be removed
fns_toberemoved <- rownames(TJcat[!TJcat$nameInTJ2level_CAs=="" & !TJcat$replaceable=="" & TJcat$v==2,])
length(fns_toberemoved)
listEntries_toberemoved <- TJcat$nameInTJ2level_CAs[!TJcat$nameInTJ2level_CAs=="" & !TJcat$replaceable=="" & TJcat$v==2]
length(listEntries_toberemoved)

TJ2level_CAs[listEntries_toberemoved] <- NULL
TJcat[fns_toberemoved,]$nameInTJ2level_CAs <- ""

## check whether it still works
pos <- 0
for (i in 1:nrow(TJcat)){
  if (TJcat$v[i]==3) next
  if (TJcat$nameInTJ2level_CAs[i]=="" & TJcat$code[i]=="") next
  pos <- pos + 1
  print(pos)
  print(rownames(TJcat)[i])
  print(TJcat$N[i]>=nrow(tjCA(TJcat$t[i], TJcat$k[i], 2)))
}

## removed files from TJ2level_CAs (August 15 2025)
## some files were still there though they names were NA
TJ2level_CAs[which(is.na(names(TJ2level_CAs)))] <- NULL
# file.remove(paste0("D:/rtests/CAs/inst/extdata/TJ2level/", fns_toberemoved))
save(TJcat, file="D:/rtests/CAs/data/TJcat.rda", compress="xz")
save(TJ2level_CAs, file="D:/rtests/CAs/data/TJ2level_CAs.rda", compress="xz")

## check whether existing constructions need replacement
settings <- TJcat[which(!TJcat$replaceable==""),]
settings$replaceable2 <- ""
for (i in 1:nrow(settings)){
  hilf <- bestN(settings$t[i], settings$k[i], 2)
  if (hilf<=settings$N[i]) settings$replaceable2[i] <- names(hilf)
}
table(settings$replaceable, settings$replaceable2)
## CS_CK and PALEY exchangeable, rest perfectly aligned
## perfectly aligned!
