## original code for reading the 2-level files in Feb 2025

# https://www.tamps.cinvestav.mx/~oc/CA/v2/t3/N54k968v2%5E968t3.ca

## t and v changed manually (t=2 with v=3 (for CAEX), t=3,4,5,6 for v=2 (for TJ2level)
t <- 4; v <- 2
tabs <- read.table(paste0("TJ_t", t, "v", v, ".txt"), head=TRUE)[[1]]
tabs <- gsub("CA(", "", tabs, fixed=TRUE)
tabs <- gsub(")≤", ",", tabs, fixed=TRUE)
(tabs <- t(sapply(strsplit(tabs, split=",", fixed=TRUE), function(obj) as.numeric(obj))))
colnames(tabs) <- c("t", "k", "v", "N")

for (i in 1:nrow(tabs)){
  k <- tabs[i,2]; N <- tabs[i,4]
  hilf <- readLines(con=paste0("https://www.tamps.cinvestav.mx/~oc/CA/v",v,"/t",t,"/N",
                               N, "k", k, "v", v, "^", k, "t", t, ".ca"))
  writeLines(hilf, con=paste0("ca.",t,".",v,"^",k,".txt"))
}
