pfad_extdata <- system.file("extdata", package="CAs")
t2 <- list.files(path=paste(pfad_extdata, "CKRS/t=2", sep="/"))
t3 <- list.files(path=paste(pfad_extdata, "CKRS/t=3", sep="/"))
t4 <- list.files(path=paste(pfad_extdata, "CKRS/t=4", sep="/"))
t5 <- list.files(path=paste(pfad_extdata, "CKRS/t=5", sep="/"))
t6 <- list.files(path=paste(pfad_extdata, "CKRS/t=6", sep="/"))


CKRScat <- rbind(data.frame(t=2, fn=t2), data.frame(t=3, fn=t3), data.frame(t=4, fn=t4),
                   data.frame(t=5, fn=t5), data.frame(t=6, fn=t6))
CKRScat[,"fn"] <- sapply(strsplit(CKRScat[,"fn"], ".", fixed=TRUE), function(obj) obj[[1]])

hilf <- CKRScat[,"fn"]
CKRScat <- CKRScat[-union(grep("_B", hilf), grep("_C", hilf)),]

crosssum <- grep("x", CKRScat$fn, fixed=TRUE)

other <- setdiff(1:nrow(CKRScat), crosssum)

CKRScat[other,]

table(colbournBigFrame[grep("CKRS", colbournBigFrame$Source),]$Source)
CKRScat$N <- sapply(CKRScat$fn, function(obj){
  hilf <- substr(obj, 4, nchar(obj))
  N <- strsplit(hilf, ";", fixed=TRUE)[[1]][1]
  print(N)
  if (length(grep("x", N))>0) N <- prod(as.numeric(unlist(strsplit(N, "x", fixed=TRUE)))) else
    N <- as.numeric(N)
  N
})
CKRScat$method <- sapply(CKRScat$fn, function(obj){
  hilf <- substr(obj, 4, nchar(obj))
  N <- strsplit(hilf, ";", fixed=TRUE)[[1]][1]
  if (length(grep("x", N))>0) method<-"crosssum" else method <- "other"
  method
})
CKRScat$v <- sapply(CKRScat$fn, function(obj){
  hilf <- substr(obj, 4, nchar(obj)-1)
  as.numeric(strsplit(hilf, ",", fixed=TRUE)[[1]][3])
})

CKRScat$k <- sapply(CKRScat$fn, function(obj){
  hilf <- substr(obj, 4, nchar(obj)-1)
  as.numeric(strsplit(hilf, ",", fixed=TRUE)[[1]][2])
})

CKRScat$Source <- sapply(1:nrow(CKRScat),
                              function(obj){
    pick <- which(colbournBigFrame$t==CKRScat[obj,]$t &
                   colbournBigFrame$N==CKRScat[obj,]$N &
                   colbournBigFrame$v==CKRScat[obj,]$v &
                   colbournBigFrame$k==CKRScat[obj,]$k)
    ifelse(length(pick)>0, colbournBigFrame$Source[pick], NA)
})

CKRScat <- CKRScat[,c("t","v","k","N","fn")]
rownames(CKRScat) <- NULL

CKRS_CAs <- lapply(1:nrow(CKRScat), function(obj){
  fn <- paste0(pfad_extdata,"/CKRS/t=", CKRScat$t[obj], "/", CKRScat$fn[obj], ".txt")
  readCA(fn, ninstruct=0, origin="CKRS", nosep = TRUE)
})
## warnings about incomplete last lines can be ignored
names(CKRS_CAs) <- CKRScat$fn

## dimensions are OK
hilf <- t(sapply(CKRS_CAs, dim))
all(hilf[,1]==CKRScat$N)
all(hilf[,2]==CKRScat$k)

## check for NAs, none found
which(sapply(CKRS_CAs, function(obj) sum(is.na(obj))>0))

## check the coverage
for (i in 1:nrow(CKRScat)){
   print(i)
  if (!all(coverage(CKRS_CAs[[i]], CKRScat$t[i])==1)) break
}
## no complaints, last one manually checked:
coverage(CKRS_CAs[[nrow(CKRScat)]], CKRScat$t[nrow(CKRScat)])

## is it possible to bring these to PCA structure ?
