
# first level: t; second level: v; third level: N
## lineages for the CA EXtender construction of Torres-Jiménez et al. (2021)
CAEX_lineages <- list(
  '2'=list(
    '3'=list(
      '11'="ca11.2.3.5",
      '12'="ca12.2.3.7",
      '13'="ca13.2.3.9",
      '14'="ca14.2.3.10",
      '15'="ca15.2.3.20",
      '16'="ca16.2.3.21",
      '17'="ca17.2.3.29",
      '18'="ca18.2.3.46",
      '19'="ca19.2.3.49",
      '20'="ca20.2.3.63",
      '21'="ca21.2.3.93",
      '22'="ca22.2.3.107",
      '23'="ca23.2.3.138",
      '24'="ca24.2.3.199",
      '25'="ca25.2.3.216",
      '26'="ca26.2.3.288",
      '27'="ca27.2.3.435",
      '28'="ca28.2.3.449",
      '29'="ca29.2.3.610",
      '30'="ca30.2.3.878",
      '31'="ca31.2.3.964",
      '32'="ca32.2.3.1308",
      '33'="ca33.2.3.1964",
      '34'=list(method="CAxCA", one="ca18.2.3.46", two="ca18.2.3.46"),
      '35'="ca35.2.3.2700",
      '36'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca21.2.3.93"),
      '37'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca22.2.3.107"),
      '38'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca23.2.3.138"),
      '39'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca24.2.3.199"),
      '40'=list(method="PCAxPCA", one="ca21.2.3.93", two="ca22.2.3.107"),
      '41'=list(method="PCAxPCA", one="ca21.2.3.93", two="ca23.2.3.138"),
      '42'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca27.2.3.435"),
      '43'=list(method="CAxCA", one="ca18.2.3.46", two="ca27.2.3.435"),
      '44'=list(method="PCAxPCA", one="ca23.2.3.138", two="ca24.2.3.199"),
      '45'=list(method="PCAxPCA", one="ca21.2.3.93", two="ca27.2.3.435"),
      '46'=list(method="PCAxPCA", one="ca22.2.3.107", two="ca27.2.3.435"),
      '47'=list(method="PCAxPCA", one="ca23.2.3.138", two="ca27.2.3.435"),
      '48'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca33.2.3.1964"),
      '49'=list(method="CAxCA", one="ca18.2.3.46", two="ca33.2.3.1964"),
      '50'=list(method="PCAxPCA", one="ca18.2.3.46", two="ca35.2.3.2628")
    )))

## reading the files that were downloaded from Torres-Jiménez website on Feb 6, 2025
toberead <- unlist(CAEX_lineages[[1]][[1]][which(sapply(CAEX_lineages[[1]][[1]],
                                                        is.character))])
fns <- strsplit(toberead, ".", fixed=TRUE)
fns <- sapply(fns, function(obj) paste0("CA.2.3^", obj[4], ".txt"))
## temporary location, move elsewhere!
#pfad <- "D:\\Users\\groemping\\ownCloud\\Documents\\Papers_DoE\\IWCT14"
pfad <- paste(system.file("extdata", package="CAs"), "CAEX", sep="/")
## the actual arrays
CAEX_CAs <- lapply(fns, function(obj) readCA(paste(pfad, obj, sep="/"),
                                 flexible.symbols = 3))
names(CAEX_CAs) <- toberead
CAEX_CAs <- lapply(CAEX_CAs, function(obj){
  class(obj) <- c("ca", class(obj))
  attr(obj, "origin") <- "Torres-Jimenez et al. 2021, CA EXtender"
  obj
})

# dimsAchieved <- sapply(11:50, function(obj) dim(CAEX(N=obj)))
# cbind(dimsAchieved, TJcat[which(TJcat[,"t"]==2)[-1],]))
#                             t      k v  N
# ca.2.3^5.txt      11      5 2      5 3 11
# ca.2.3^7.txt      12      7 2      7 3 12
# ca.2.3^9.txt      13      9 2      9 3 13
# ca.2.3^10.txt     14     10 2     10 3 14
# ca.2.3^20.txt     15     20 2     20 3 15
# ca.2.3^21.txt     16     21 2     21 3 16
# ca.2.3^29.txt     17     29 2     29 3 17
# ca.2.3^46.txt     18     46 2     46 3 18
# ca.2.3^49.txt     19     49 2     49 3 19
# ca.2.3^63.txt     20     63 2     63 3 20
# ca.2.3^93.txt     21     93 2     93 3 21
# ca.2.3^107.txt    22    107 2    107 3 22
# ca.2.3^138.txt    23    138 2    138 3 23
# ca.2.3^199.txt    24    199 2    199 3 24
# ca.2.3^216.txt    25    216 2    216 3 25
# ca.2.3^288.txt    26    288 2    288 3 26
# ca.2.3^435.txt    27    435 2    435 3 27
# ca.2.3^449.txt    28    449 2    449 3 28
# ca.2.3^610.txt    29    610 2    610 3 29
# ca.2.3^878.txt    30    878 2    878 3 30
# ca.2.3^964.txt    31    964 2    964 3 31
# ca.2.3^1308.txt   32   1308 2   1308 3 32
# ca.2.3^1964.txt   33   1964 2   1964 3 33
# ca.2.3^2116.txt   35   2116 2   2116 3 34  *** productCA1
# ca.2.3^2700.txt   35   2700 2   2700 3 35
# ca.2.3^3918.txt   36   3918 2   3918 3 36
# ca.2.3^4457.txt   37   4457 2   4457 3 37
# ca.2.3^5763.txt   38   5763 2   5763 3 38
# ca.2.3^8329.txt   39   8329 2   8329 3 39
# ca.2.3^9207.txt   40   9207 2   9207 3 40
# ca.2.3^11898.txt  41  11898 2  11898 3 41
# ca.2.3^17895.txt  42  17895 2  17895 3 42
# ca.2.3^20010.txt  44  20010 2  20010 3 43  *** productCA1
# ca.2.3^25317.txt  44  25317 2  25317 3 44
# ca.2.3^37071.txt  45  37071 2  37071 3 45
# ca.2.3^42174.txt  46  42174 2  42174 3 46
# ca.2.3^54531.txt  47  54531 2  54531 3 47
# ca.2.3^80789.txt  48  80789 2  80789 3 48
# ca.2.3^90344.txt  50  90344 2  90344 3 49  *** productCA1
# ca.2.3^112770.txt 50 112770 2 112770 3 50
