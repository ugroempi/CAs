## copy-pasted from Lobb et al. 2012 (LCDST)
##   replaced oo_0 with 100+v-1, ..., oo_f-1 with 100+v-f
## (initially?) only included cases that are used with group Z_v-f
## first list level: v
## second list level: k

LCDSTStarters <- list(
  '3' = list('5'=c(102,0,1,1,0)),
  '4' = list('6'=c(103,0,1,1,2,1),
             '7'=c(103,2,1,2,2,1,2),
             '9'=c(103,102,1,0,1,1,102,103,0)),
  '5' = list('7'=c(104,2,3,1,1,0,1),
             # '8'=c(104,2,3,1,1,0,1),  ## distinct, k-1 elements, same as seven
             # can omit, as covered in CS_MS
             '10'=c(103,103,2,1,104,1,104,0,1,1)),
  '6' = list('9'=c(105,2,4,4,2,1,2,1,3)),
  '7' = list('10'=c(106,0,1,0,3,5,4,2,0,0),
             '13'=c(106,3,106,105,3,2,1,105,0,2,0,1,1),
             '14'=c(106,1,4,0,2,1,4,4,4,4,106,105,0,105),
             '15'=c(105,105,3,106,4,1,2,106,2,2,NA,2,2,1,4),
             '17'=c(104,104,1,105,3,1,2,106,1,106,3,106,3,2,2,105,3)),
  '8' = list('11'=c(107,0,3,0,6,4,1,3,0,1,1),
             ## 13 has to be figured out
             ## is there are group for 8 that does not have a single primitive only?
             ## x x^2 x xy xy 107 1 x^y x^2 107 106 1 106
             ##           in {x, yx^3 = y^2 = 1, yxy = x^2}
             '15'=c(0,1,3,3,1,2,5,106,0,0,5,106,107,0,107),
             '16'=c(107,106,2,107,1,106,5,5,1,1,2,0,2,5,3,2),
             '18'=c(107,107,1,3,3,106,1,106,0,4,0,105,1,0,105,0,2,0),
             '19'=c(107,107,0,3,106,0,105,4,4,0,NA,0,2,1,105,0,106,3,4),
             '20'=c(107,107,1,106,0,0,1,106,0,105,4,3,1,3,105,3,0,2,3,0)),
  '9' = list('13'=c(108,7,0,3,2,4,1,4,4,0,1,1,7),
             '16'=c(3,108,3,2,107,6,5,5,1,6,3,5,6,108,4,107),
             '17'=c(2,4,108,3,4,6,6,3,6,6,4,0,108,6,5,107,107),
             '19'=c(108,108,3,107,2,106,3,2,2,107,3,0,1,5,3,1,3,106,0),
             '20'=c(108,108,5,4,107,2,107,1,2,106,5,1,NA,1,5,2,106,3,3,5),
             '21'=c(108,108,3,107,4,2,3,107,3,106,1,3,2,5,4,4,NA,1,0,106,2),
             '22'=c(108,108,4,107,3,2,2,107,3,106,5,1,2,5,106,5,NA,0,3,3,4,2)),
  '10' = list( ## omit 14 (group Z3xZ3), 16 (distinct, k-1, not better than MS), 25 to 27 (f=4)
             '15'=c(109,4,6,1,6,3,7,8,2,0,1,5,4,6,6),
             '18'=c(109,5,3,3,109,3,6,3,4,4,0,3,5,0,7,108,2,108),
             '20'=c(109,109,5,108,0,5,4,108,5,107,0,3,5,2,107,0,1,4,4,4),
             '21'=c(109,108,3,107,0,6,4,107,0,108,109,6,5,3,0,2,5,0,0,1,0),
             '22'=c(109,109,3,108,4,0,6,108,2,107,3,2,3,5,107,4,NA,5,3,3,3,0),
             '23'=c(109,109,4,108,6,6,3,1,4,108,5,4,107,4,NA,4,107,5,0,1,6,NA,4),
             '25'=c(109,109,4,108,1,107,3,1,2,108,5,0,2,106,2,1,107,2,106,2,5,106,5,3,3),
             '26'=c(109,0,1,108,109,4,107,3,5,5,2,0,106,107,106,106,1,5,4,5,1,0,2,1,108,5),
             '27'=c(109,109,4,108,1,4,3,108,5,107,4,2,106,0,1,5,5,5,1,107,1,2,0,3,0,106,1)),
  '11' = list('16'=c(110,2,6,7,2,1,4,5,3,3,6,2,9,9,0,2),
              '17'=c(110,1,0,5,5,8,2,4,4,1,2,6,2,5,3,2,0),
              '19'=c(3,2,7,110,5,4,8,0,0,109,110,4,7,1,3,0,3,1,109),
              '20'=c(110,109,3,110,1,109,6,3,5,2,5,4,4,8,6,2,5,1,2,2),
              '22'=c(110,110,5,109,4,7,1,109,3,108,0,1,0,5,108,1,3,7,4,2,7,7),
              '23'=c(110,109,1,108,2,6,6,108,7,109,110,3,2,1,2,6,4,1,0,2,5,0,4),
              '25'=c(110,110,0,109,0,108,0,3,1,109,2,4,4,107,1,2,108,6,107,5,0,107,5,2,1),
              '26'=c(110,1,5,109,110,3,108,4,0,0,5,0,107,108,107,107,1,1,5,3,2,1,0,1,109,0),
              '27'=c(110,110,0,109,4,3,4,109,4,108,5,4,107,2,0,2,3,6,6,108,4,3,3,0,5,107,6))
)

# initial implementation for f<=3 could work with eCAN
eCANhilf <- function(k,v){
  if (v==1) return(1) else return(eCAN(2,k,v)$CAN)
}
eCANvec <- Vectorize(eCANhilf)
bestNvec <- function(k, v, ...) mapply(function(obj1, obj2)
    unname(ifelse(obj2==1, 1, bestN(2,obj1,obj2, ...))),
    k, v)


## obtain realistic N via applying function bestN for the CA(2,k,f)
## implies that updates are needed if best designs for relevant cases are added
##     to the package
LCDSTCombis <- data.frame(
  v=(v <- rep(as.numeric(names(LCDSTStarters)), lengths(LCDSTStarters))),
  k=(k <- unlist(lapply(LCDSTStarters, function(obj) as.numeric(names(obj))))),
  f=(f <- unlist(lapply(LCDSTStarters, function(obj) sapply(obj,
                        function(obj2) length(unique(obj2[which(obj2>100)])))))),
  NAs=(NAs <- unlist(lapply(LCDSTStarters, function(obj) sapply(obj,
                                              function(obj2) sum(is.na(obj2)))))),
  N = (v-f)*k + bestNvec(k,f),
  N_withoutInternet = (v-f)*k + bestNvec(k,f, internet=FALSE),
  N_eCANbased = (v-f)*k + eCANvec(k,f)
)

## verify coverage and dimensions
for (i in 1:nrow(LCDSTCombis)){
  print(i)
  if (!all(coverage(D<-CS_LCDST(LCDSTCombis$k[i], LCDSTCombis$v[i]), 2)==1)) break
  if (!(nrow(D)==LCDSTCombis$N[i] && ncol(D)==LCDSTCombis$k[i])) break
  if (f[i]>3) print(paste0("excess size: ", LCDSTCombis$N[i] - LCDSTCombis$N_eCANbased[i]))
}

