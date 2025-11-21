## the string is copy-pasted from
## https://github.com/robbywalker/ca-phf-research/blob/master/results/results.txt

## the actual arrays sometimes have a few more columns
## than in the strings below

string <-
  "CAN(3,22,3) <= 75
CAN(3,33,3) <= 99
CAN(3,52,3) <= 123
CAN(3,120,3) <= 171

CAN(3,16,4) <= 124
CAN(3,34,4) <= 184
CAN(3,60,4) <= 244
CAN(3,110,4) <= 304
CAN(3,160,4) <= 364

CAN(3,44,5) <= 365
CAN(3,93,5) <= 485
CAN(3,160,5) <= 605

CAN(3,32,7) <= 679
CAN(3,74,7) <= 1015
CAN(3,150,7) <= 1351

CAN(3,30,8) <= 1016
CAN(3,91,8) <= 1520
CAN(3,200,8) <= 2024

CAN(3,41,9) <= 1449
CAN(3,113,9) <= 2169
CAN(3,225,9) <= 2889

CAN(3,49,11) <= 2651
CAN(3,146,11) <= 3971

CAN(3,58,13) <= 4381
CAN(3,200,13) <= 6565

CAN(4,9,2) <= 44
CAN(4,12,2) <= 58

CAN(4,10,3) <= 159
CAN(4,22,3) <= 315
CAN(4,28,3) <= 393
CAN(4,35,3) <= 471
CAN(4,46,3) <= 549

CAN(4,13,4) <= 508
CAN(4,20,4) <= 760
CAN(4,30,4) <= 1012
CAN(4,42,4) <= 1264

CAN(4,14,5) <= 1245
CAN(4,23,5) <= 1865
CAN(4,35,5) <= 2485
CAN(4,62,5) <= 3105

CAN(4,18,7) <= 4795
CAN(4,30,7) <= 7189
CAN(4,54,7) <= 9583

CAN(4,19,8) <= 8184
CAN(4,36,8) <= 12272

CAN(4,21,9) <= 13113
CAN(4,39,9) <= 19665

CAN(4,23,11) <= 29271

CAN(5,10,2) <= 62
CAN(5,12,2) <= 92
CAN(5,14,2) <= 122
CAN(5,16,2) <= 152
CAN(5,17,2) <= 182

CAN(5,10,3) <= 483
CAN(5,13,3) <= 723
CAN(5,16,3) <= 963
CAN(5,19,3) <= 1203

CAN(5,11,4) <= 2044
CAN(5,15,4) <= 3064

CAN(6,10,3) <= 1455

CAN(7,9,2) <= 254
CAN(7,11,2) <= 380
CAN(7,12,2) <= 506
CAN(7,13,2) <= 758"

## figure out the number of improvements obtainable
zeilen <- strsplit(string, "\n", fixed=TRUE)[[1]]
zeilen <- zeilen[nchar(zeilen)>0]
length(zeilen)
code <- sapply(strsplit(zeilen, " "), function(obj) obj[[1]])
code <- gsub("CAN", "bestN", code)
jetzt <- sapply(code[1:60], function(obj) eval(parse(text=obj)))
neu <- as.numeric(sapply(strsplit(zeilen, " "), function(obj) obj[[3]]))[1:60]
length(besser <- which(neu<jetzt))
sort(jetzt[besser]-neu[besser])
## 36 improvements, from a few runs up to several thousand runs better
code[besser]
## namenbesser was manually populated,
## looking up the cphf names from the repository
namenbesser <- c(
  "4-37-3-3.cphf", #bestN(3,33,3)"
  "4-64-4-3.cphf", # bestN(3,60,4)
  "5-120-4-3.cphf", # bestN(3,110,4)
  "3-48-5-3.cphf", #bestN(3,44,5)",
  "4-95-5-3.cphf", # 93
  "2-32-7-3.cphf",
  "3-81-7-3.cphf", #bestN(3,74,7)",
  "3-91-8-3.cphf",
  "2-41-9-3.cphf",  "3-113-9-3.cphf",
  "2-50-11-3.cphf", "3-146-11-3.cphf",
  "2-59-13-3.cphf", "3-200-13-3.cphf",
  "2-13-4-4.cphf",  "3-20-4-4.cphf",
  "4-31-4-4.cphf", # bestN(4,30,4)",
  "5-42-4-4.cphf",
  "2-15-5-4.cphf", #bestN(4,14,5)",
  "3-24-5-4.cphf", # bestN(4,23,5)"
  "4-37-5-4.cphf", # bestN(4,35,5)",
  "5-62-5-4.cphf", "2-18-7-4.cphf",
  "3-30-7-4.cphf", "4-54-7-4.cphf",
  "2-20-8-4.cphf", # "bestN(4,19,8)",
  "3-36-8-4.cphf",  "2-21-9-4.cphf",
  "3-39-9-4.cphf",  "2-23-11-4.cphf",
  "2-10-3-5.cphf",  "3-13-3-5.cphf",
  "4-16-3-5.cphf",  "2-11-4-5.cphf",
  "3-15-4-5.cphf",  "2-10-3-6.cphf")

pfadGithub <- "https://raw.githubusercontent.com/robbywalker/ca-phf-research/master/results"
## pfad <- paste0(pfadGithub, "/", nam)
## cphf <- read.csv(pfad, skip=1, header=FALSE)

# ## download test
# download.file(url=paste0(pfadGithub, "/", namenbesser[1]),
#              destfile = paste0("D:/rtests/CAsbak/arrays/WC/", namenbesser[1]))
#
 # for (fn in namenbesser){
 #   download.file(url=paste0(pfadGithub, "/", fn),
 #                 destfile = paste0("D:/rtests/CAsbak/arrays/WC/", fn))
 # }

transform <- function(cphf_string, s){
  ## assumes that the symbol q-1 occurs, which should always be the case, I hope
  ## s is the length of the tuples
  ## CPHF columns are denoted with numbers
  ## 0, ..., q^s-1
  hilf <- strsplit(cphf_string, " ")
  hilf <- lapply(hilf, function(obj)
    sapply(obj, function(obj2)
      unlist(strsplit(obj2, ""))))
  alphabet <- sort(unique(unlist(hilf)))
  if ("A" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="A")] <- "10"; obj})
  if ("B" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="B")] <- "11"; obj})
  if ("C" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="C")] <- "12"; obj})
  if ("D" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="D")] <- "13"; obj})
  if ("E" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="E")] <- "14"; obj})
  if ("F" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="F")] <- "15"; obj})
  alphabet <- sort(as.numeric(unique(unlist(hilf))))
  q <- max(alphabet) + 1
  t(sapply(hilf, function(obj){
    dd <- dim(obj)
    obj <- matrix(as.numeric(obj), nrow=dd[1])
    t(q^((s - 1) : 0)) %*% obj
  }))
}

transformCPHF <- function(cphf_string, s){
  ## assumes that the symbol q-1 occurs, which should always be the case, I hope
  ## s is the length of the tuples after rbinding them with 1
  ##          (i.e., before, they are from the Walker/Colbourn paper
  ##           or the github repo)
  ## CPHF columns are denoted with numbers
  ## 0, ..., q^s-1
  hilf <- strsplit(cphf_string, " ")
  hilf <- lapply(hilf, function(obj)
    sapply(obj, function(obj2)
      unlist(strsplit(obj2, ""))))
  alphabet <- sort(unique(unlist(hilf)))
  if ("A" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="A")] <- "10"; obj})
  if ("B" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="B")] <- "11"; obj})
  if ("C" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="C")] <- "12"; obj})
  if ("D" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="D")] <- "13"; obj})
  if ("E" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="E")] <- "14"; obj})
  if ("F" %in% alphabet) hilf <- lapply(hilf,
                                        function(obj) {obj[which(obj=="F")] <- "15"; obj})
  alphabet <- sort(as.numeric(unique(unlist(hilf))))
  q <- max(alphabet) + 1
  t(sapply(hilf, function(obj){
    dd <- dim(obj)
    obj <- rbind(matrix(as.numeric(obj), nrow=dd[1]),1)
    t(q^((s - 1) : 0)) %*% obj
  }))
}

## read from github repo
## and subsequently transformed
## the Galois fields are like for SMC, except for q=9,
##    where the Galois field must be taken from
##    https://github.com/robbywalker/ca-phf-research/blob/master/tabu/math/galois_field.cpp
WC_SCPHFs <- list(
   '3'=vector(mode="list", 0),
   '4'=vector(mode="list", 0),
   '5'=vector(mode="list", 0),
   '6'=vector(mode="list", 0)
)
for (fn in namenbesser){
  cphf <- read.csv(paste0("D:/rtests/CAsbak/arrays/WC/", fn),
                   skip=1, header=FALSE)
  params <- read.table(paste0("D:/rtests/CAsbak/arrays/WC/", fn),
                       nrows=1, header=FALSE)
  names(params) <- c("ncphf", "k", "v", "t")

  (cphf <- transform(unlist(cphf), params$t-1))
  # ## coverage OK
  # if (!all(coverage(SMC(cphf, params$v, params$t), params$t, parallel=6)==1)){
  #   message("error for ", fn)
  #   next
  # }
  WC_SCPHFs[[as.character(params$t)]][[as.character(params$v)]][[as.character(params$k)]] <-
    cphf
}

for (t in 3:6){
#for (t in 4:4){
  tliste <- WC_SCPHFs[[as.character(t)]]
  print(paste0("t=", t))
  # for (v in names(tliste)){
  for (v in "9"){
    vliste <- tliste[[v]]
    print(paste0("v=", v))
    for (k in names(vliste)){
      print(k)
      hilf <- SMC(vliste[[k]], as.numeric(v), t, type="2018")
      erg <- coverage(hilf, t, parallel=6)
      if (all(erg==1)) print("ok") else print(erg)
    }
  }
}

WCcat <- data.frame(t=numeric(0), v=numeric(0), k=numeric(0), N=numeric(0))
for (t in 3:6)
  for (v in as.numeric(names(WC_SCPHFs[[as.character(t)]]))){
    for (k in as.numeric(names(WC_SCPHFs[[as.character(t)]][[as.character(v)]]))){
      kc <- as.character(k); vc <- as.character(v); tc <- as.character(t)
      WCcat <- rbind(WCcat,
      c(t=t, v=v, k=k,
        N=nrow(SMC(WC_SCPHFs[[tc]][[vc]][[kc]], v, t))))
    }
  }
colnames(WCcat) <- c("t","v", "k", "N")

save("WCcat", file="d:/rtests/CAsbak/arrays/WCcat.rda", compress="xz")
save("WC_SCPHFs", file="d:/rtests/CAs/data/WC_SCPHFs.rda", compress="xz")

SCPHFcat <- rbind(SCPHFcat, cbind(WCcat, type="2009"))
SCPHFcat <- SCPHFcat[ord(SCPHFcat),]
rownames(SCPHFcat) <- NULL
save("SCPHFcat", file="d:/rtests/CAs/data/SCPHFcat.rda", compress="xz")
