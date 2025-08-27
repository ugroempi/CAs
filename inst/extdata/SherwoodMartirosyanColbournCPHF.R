
transform <- function(cphf_string, s){
  ## assumes that the symbol q-1 occurs
  ## pow is the power so that the number of levels is v=q^s
  ## as 0, ..., v-1
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
  ## assumes that the symbol q-1 occurs
  ## s is the power so that the number of levels is v=q^s
  ## as 0, ..., v-1
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

## copy-pasted from Sherwood, Martirosyan, Colbourn (2006)
  ############## strength t=3 ####################################
#CPHF(3; 20, 9, 3)  # v=3 # paper: N=75, k=20
SMC_CPHFs <- list(
   '3'=vector(mode="list", 0),
   '4'=vector(mode="list", 0)
)
cphf <- c("00 01 02 10 11 12 20 21 22 00 01 02 10 11 12 20 21 22 00 01",
          "00 01 02 10 11 22 20 21 12 22 20 21 12 10 01 02 00 11 20 00",
          "00 01 10 02 11 21 22 12 20 10 11 01 00 12 20 02 22 21 20 21")
(cphf <- transform(cphf, 2))
## coverage OK
SMC_CPHFs$`3`[["3"]] <- list("20"=cphf)

#CPHF(2; 16, 16, 3)  # v=4 # paper: N=124, k=16
cphf <- c("00 01 02 03 10 11 12 13 20 21 22 23 30 31 32 33",
          "00 01 10 11 02 03 12 13 21 20 31 30 23 22 33 32")
(cphf1 <- transform(cphf, 2))
## coverage OK

#CPHF(3; 28, 16, 3)  # v=4 # paper: N=184, k=28
cphf <- c("00 01 02 03 10 11 12 13 20 21 22 23 30 31 32 33 00 01 02 03 10 11 12 13 20 21 22 23",
          "00 01 02 03 10 11 12 13 20 21 22 23 30 31 33 32 11 10 22 23 01 00 20 21 32 33 02 03",
          "00 01 10 11 02 03 12 13 21 20 31 30 23 22 32 33 11 10 01 00 31 30 13 12 22 23 02 03")
(cphf2 <- transform(cphf, 2))
## coverage OK
SMC_CPHFs$`3`[["4"]] <- list("16"=cphf1, "28"=cphf2)

#CPHF(2; 24, 25, 3)  # v=5 # paper: N=245, k=24
cphf <- c("00 01 02 03 04 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43",
          "00 01 10 11 23 04 42 03 13 24 14 02 31 41 43 22 30 21 40 32 12 44 33 34")
(cphf <- transform(cphf, 2))
## coverage OK
SMC_CPHFs$`3`[["5"]] <- list("24"=cphf)

#CPHF(2; 31, 49, 3)  # v=7 # paper: N=679, k=31
cphf <- c("00 01 02 03 04 05 06 10 11 12 13 14 15 16 20 21 22 23 24 25 26 30 31 32 33 34 35 36 40 41 42",
          "00 01 10 11 23 25 63 02 03 12 13 20 60 65 22 43 46 62 55 16 53 50 52 61 42 15 31 40 34 35 41")
(cphf <- transform(cphf, 2))
## coverage OK
SMC_CPHFs$`3`[["7"]] <- list("31"=cphf)

#CPHF(2; 40, 64, 3)  # v=8 # paper: N=1016, k=40
cphf <- c("00 01 02 03 04 05 06 07 10 11 12 13 14 15 16 17 20 21 22 23 24 25 26 27 30 31 32 33 34 35 36 37 40 41 42 43 44 45 46 47",
          "00 01 10 11 24 25 34 35 02 03 12 13 26 27 36 37 21 20 31 30 05 04 15 14 23 22 33 32 07 06 17 16 51 50 41 40 75 74 65 64")
(cphf <- transform(cphf, 2))
## this uses a different primitive polynomial than package lhs for the Galois field:
##     lhs uses x^3 = 1 + x^2 (RHS in xton)
##         which is equivalent to x^3 - x^2 - 1 = 0
##         which is equivalent to x^3 + x^2 + 1 = 0 in terms of GF(2),
##         i.e., primitive element 13
##     presumably, the paper uses x^3 + x + 1 = 0, i.e., primitive element 11
##

## with the correct galois field, it works with all 40 columns
## with the default one, only with the first 32 columns
SMC_CPHFs$`3`[["8"]] <- list("40"=cphf)

#CPHF(2; 39, 81, 3)  # v=9 # paper: N=1449, k=39
cphf <- c("00 01 02 03 04 05 06 07 08 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 40 41 42",
          "00 01 10 11 24 26 53 57 83 02 03 12 13 20 25 46 48 66 14 15 71 64 30 80 76 85 68 74 67 08 51 41 47 78 62 05 21 52 84")
(cphf <- transform(cphf, 2))
## this also does not work with the primitive polynomial of package lhs for the Galois field:
##     lhs uses x^2 = 1 + 2x (RHS in xton)
##         which is equivalent to x^2 - 2x - 1 = 0
##         which is equivalent to x^2 + x + 2 = 0 in terms of GF(3),
##         i.e., primitive element 14
## the suitable characteristic polynomial is x^2+2x+2
##
SMC_CPHFs$`3`[["9"]] <- list("39"=cphf)


#CPHF(2; 44, 121, 3) # v=11 # paper: N=2651, k=44
cphf <- c("00 01 02 03 04 05 06 07 08 09 0A 10 11 12 13 14 15 16 17 18 19 1A 20 21 22 23 24 25 26 27 28 29 2A 30 31 32 33 34 35 36 37 38 39 3A",
          "00 01 10 11 24 28 43 49 83 89 A4 02 03 12 13 26 2A 40 45 80 85 A6 14 21 20 04 0A A7 70 77 19 51 58 50 64 35 86 96 A0 A5 34 99 68 88")
(cphf <- transform(cphf, 2))
## coverage OK
SMC_CPHFs$`3`[["11"]] <- list("44"=cphf)


#CPHF(2; 39, 169, 3) # v=13
cphf <- c(
  "00 01 02 03 04 05 06 07 08 09 0A 0B 0C 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 20 21 22 23 24 25 26 27 28 29 2A 2B 2C",
  "00 01 10 11 24 2A 45 49 57 97 A5 A9 C4 02 03 12 13 26 2C 47 4B 59 99 A7 AB C6 14 23 20 06 0A 1C 73 84 70 8C B2 B1 9A")
(cphf <- transform(cphf, 2))
# coverage(SMC(cphf,13, 3),3, parallel = 4)
# covers

SMC_CPHFs$`3`[["13"]] <- list("39"=cphf)


############## strength t=4 ################
#CPHF(2; 8, 8, 4)  # v=2 # paper: N=30, k=8
cphf <- c(
  "000 001 010 011 100 101 110 111",
  "000 001 010 100 011 110 111 101")
(cphf1 <- transform(cphf, 3))


#CPHF(4; 11, 8, 4) # v=2 # paper: N=58, k=11
cphf <- c("000 001 010 011 100 101 110 111 000 001 010",
          "000 001 010 011 100 101 111 110 101 111 001",
          "000 001 010 011 100 110 111 101 111 010 101",
          "000 001 010 100 011 111 101 110 011 100 000")
(cphf2 <- transform(cphf, 3))

SMC_CPHFs$`4`[["2"]] <- list("8"=cphf1,
                              "11"=cphf2)

#CPHF(2; 10, 27, 4)   # v=3 # paper: N=159, k=10
cphf <- c("000 001 002 010 011 100 101 110 120 121",
          "000 001 010 100 111 101 211 122 212 112")
(cphf1 <- transform(cphf, 3))

#CPHF(3; 16, 27, 4)   # v=3 # paper: N=237, k=16
cphf <- c("000 001 002 010 011 012 020 021 100 101 102 110 111 112 211 222",
          "000 001 002 010 011 100 101 110 012 020 221 112 220 022 111 122",
          "000 001 010 100 111 101 211 120 200 121 221 201 220 210 110 102")
(cphf2 <- transform(cphf, 3))
SMC_CPHFs$`4`[["3"]] <- list("10"=cphf1,
                              "16"=cphf2)


#CPHF(2; 9, 64, 4)    # v=4 # paper: N=508, k=9
cphf <- c("000 001 002 003 010 100 110 121 131",
          "000 001 010 100 112 113 121 211 222")
(cphf <- transform(cphf, 3))
SMC_CPHFs$`4`[["4"]] <- list("9"=cphf)


#CPHF(2; 11, 125, 4)  # v=5 # paper: N=1245, k=11
cphf <- c("000 001 002 003 004 010 100 110 121 212 231",
          "000 001 010 100 111 124 142 214 241 412 421")
(cphf <- transform(cphf, 3))
SMC_CPHFs$`4`[["5"]] <- list("11"=cphf)

SCPHFcat <- data.frame(t=numeric(0), v=numeric(0), k=numeric(0), N=numeric(0))
for (t in 3:4)
  for (v in as.numeric(names(SMC_CPHFs[[as.character(t)]]))){
    for (k in as.numeric(names(SMC_CPHFs[[as.character(t)]][[as.character(v)]]))){
      kc <- as.character(k); vc <- as.character(v); tc <- as.character(t)
      SCPHFcat <- rbind(SCPHFcat,
      c(t=t, v=v, k=k,
        N=nrow(SMC(SMC_CPHFs[[tc]][[vc]][[kc]], v, t))))
    }
  }
colnames(SCPHFcat) <- c("t","v", "k", "N")

SCPHFcat$type <- "2006"

save("SCPHFcat", file="d:/rtests/CAs/data/SCPHFcat.rda", compress="xz")
## caution, is later also modified by getCL_CPHFs_fromDwyerRepo.R
save("SMC_CPHFs", file="d:/rtests/CAs/data/SMC_CPHFs.rda", compress="xz")

## check that everything is correct
for (t in 3:4)
  for (v in as.numeric(names(SMC_CPHFs[[as.character(t)]]))){
    for (k in as.numeric(names(SMC_CPHFs[[as.character(t)]][[as.character(v)]]))){
      print(paste0("t=",t, ", v=", v, ", k=", k))
      kc <- as.character(k); vc <- as.character(v); tc <- as.character(t)
      cphf <- SMC_CPHFs[[tc]][[vc]][[kc]]
      print(coverage(SMC(cphf, v, t), t, parallel=4))
          }
  }
