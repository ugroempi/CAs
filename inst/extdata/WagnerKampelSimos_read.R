## read designs from Wagner Kampel Simos 2021
## the csv files are not included in extdata,
## as they are relatively large
## and they are available at

pathSIPO <- "https://srd.sba-research.org/data/sipo/"

## file names within that small repo

fns <- c("CA(440-6,30,2).csv",
"CA(452-6,31,2).csv", "CA(461-6,32,2).csv", "CA(470-6,33,2).csv",
"CA(479-6,34,2).csv", "CA(487-6,35,2).csv", "CA(495-6,36,2).csv",
"CA(506-6,37,2).csv", "CA(512-6,38,2).csv", "CA(518-6,39,2).csv",
"CA(527-6,40,2).csv", "CA(534-6,41,2).csv", "CA(541-6,42,2).csv",
"CA(553-6,43,2).csv", "CA(561-6,44,2).csv", "CA(566-6,45,2).csv",
"CA(573-6,46,2).csv", "CA(581-6,47,2).csv", "CA(586-6,48,2).csv",
"CA(592-6,49,2).csv", "CA(599-6,50,2).csv", "CA(605-6,51,2).csv",
"CA(612-6,52,2).csv", "CA(615-6,53,2).csv", "CA(623-6,54,2).csv",
"CA(628-6,55,2).csv", "CA(634-6,56,2).csv", "CA(640-6,57,2).csv",
"CA(644-6,58,2).csv", "CA(651-6,59,2).csv", "CA(660-6,60,2).csv",
"CA(666-6,61,2).csv", "CA(672-6,62,2).csv", "CA(677-6,63,2).csv",
"CA(680-6,64,2).csv", "CA(685-6,65,2).csv", "CA(688-6,66,2).csv",
"CA(695-6,67,2).csv", "CA(699-6,68,2).csv", "CA(705-6,69,2).csv",
"CA(709-6,70,2).csv", "CA(714-6,71,2).csv", "CA(718-6,72,2).csv")

WKS_CAs <- lapply(fns, function(obj){
  readCA(paste0(pathSIPO, obj), ninstruct = 0, sep=",",
         header=FALSE, origin="Wagner, Kampel, Simos (2021) (WKS)")
  }
  )

## names are k, as this is consecutive
names(WKS_CAs) <- 30:72

fns <- sapply(strsplit(fns, "(", fixed=TRUE), function(obj) obj[2])
fns <- sapply(strsplit(fns, ")", fixed=TRUE), function(obj) obj[1])
N <- as.numeric(sapply(strsplit(fns, "-", fixed=TRUE), function(obj) obj[1]))
tkv <- sapply(strsplit(fns, "-", fixed=TRUE), function(obj) obj[2])
tkv <- matrix(as.numeric(do.call(rbind, strsplit(tkv, ","))),
              ncol=3, dimnames=list(NULL, c("t","k","v")))

WKScat <- cbind(as.data.frame(tkv), N=N)


## check that everything is OK
cbind(WKScat, t(sapply(WKS_CAs, dim)))
