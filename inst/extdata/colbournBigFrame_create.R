## https://www.public.asu.edu/~ccolbou/src/tabby/2-5-ca.html
## t=2, v=5
## the source is no longer accessible, i.e., the code would not run now

cases <- expand.grid(t=2:6, v=2:25)
colbournCatalogue <- array(data=list(NA), dim=c(5,24), dimnames=list(t=2:6, v=2:25))
for (i in 1:nrow(cases))
{ t <- cases$t[i]
  v <- cases$v[i]
  gelesen <- readLines(paste0(
     "https://www.public.asu.edu/~ccolbou/src/tabby/",t,"-",v,"-ca.html"))
  starts <- substr(gelesen,1,7)
  von <- which(starts=="<table ") + 1 ## header row
  bis <- which(starts=="</table") - 1  ## last row
  gelesen <- gelesen[von:bis]
  ncs <- nchar(gelesen)
  gelesen <- substr(gelesen, 9, ncs-10)
  gelesen <- gsub("</td><td>", ";", gelesen)
  gelesen <- read.table(text=gelesen, header=TRUE, sep=";")
  ## v=2 and t=2 has two columns only, source in general text
  if (v==2 && t==2) gelesen <- cbind(gelesen, Source="Kleitman and Spencer, or Katona")
  gelesen <- cbind(t, v, gelesen)
  colnames(gelesen) <- c("t", "v", "k", "N", "Source")
  colbournCatalogue[t-1, v-1] <- list(gelesen)
}

## combine all files into one long data frame
colbournBigFrame <- do.call(rbind, c(colbournCatalogue))