## DwyerCAsDatabase at
## https://github.com/aadwyer/CA_Database/blob/main/Repository/CA/
## read the arrays that are from personal communication with Colbourn
## (and not from Kokkala et al. 2018)
## as of July 2, 2025

pfadGithub <- "https://raw.githubusercontent.com/aadwyer/CA_Database/main/Repository/CA"

infostrings <- readLines(paste0(pfadGithub,"/_index.txt"))
infostrings <- infostrings[setdiff(1:length(infostrings), grep("Kokkala", infostrings))]
table(duplicated(infostrings))
infostrings <- unique(infostrings)
length(infostrings)   # 307

## file names
fns <- unname(sapply(infostrings, function(obj)
  strsplit(obj, " ", fixed=TRUE)[[1]][[1]]))
infos <- unname(sapply(infostrings, function(obj)
  strsplit(obj, ", construction is ", fixed=TRUE)[[1]][[2]]))

## load from Dwyer repo
## as far as they exist (as I could not find a way to list the available
##                       arrays without loading them)
pfadGithub <- "https://raw.githubusercontent.com/aadwyer/CA_Database/main/Repository/CA"
status <- logical(length=length(fns))
## this loop takes a while
for (i in 1:length(fns)){
  nam <- fns[i]
  pfad <- paste0(pfadGithub, "/", nam)
  print(pfad)
  aus <- try(readCA(pfad, ninstruct = 0, ignore.chars=c("[","]"),
              header=FALSE, sep=",",
              origin=paste0("Dwyer Github repository, ", nam, ", ",
                            infos[i])))
if ("try-error" %in% class(aus)) status[i] <- FALSE else status[i] <- TRUE
}

table(status)

## restrict to available arrays
fns <- fns[status]
infos <- infos[status]

## omit CA_4998_6_7_4.txt
## as it is worse than the OA for the same setting

omit <- which(fns=="CA_4998_6_7_4.txt")
fns <- fns[-omit]
infos <- infos[-omit]

## create a data frame with all information
hilf <- do.call(rbind, strsplit(fns, "_", fixed=TRUE))
DWYERcat <- data.frame(t=as.numeric(hilf[,3]),
                     v=as.numeric(sapply(strsplit(hilf[,5], ".", fixed=TRUE), function(obj) obj[1])),
                     k=as.numeric(hilf[,4]),
                     N=as.numeric(hilf[,2]),
                     Source=infos)
row.names(DWYERcat) <- fns
DWYERcat <- DWYERcat[DoE.base::ord(DWYERcat),]

## test that it works
dwyerCA(4, 7, 5)
