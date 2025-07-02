## code for reading the 2-level files available in the extdata folder

## Author: Ulrike Groemping

# original location of the files in Feb 2025
# https://www.tamps.cinvestav.mx/~oc/CA/v2/t3/N54k968v2%5E968t3.ca
# source of the files: Jose Torres-Jimenez homepage

pfad_extdata <- system.file("extdata", package="CAs")
fileloc <- paste(pfad_extdata, "TJ2level", sep="/")
toberead <- list.files(fileloc)
fns <- rownames(TJcat)[which(TJcat[,"v"]==2 & !TJcat$nameInTJ2level_CAs=="")]
length(setdiff(toberead, fns))

TJ2level_CAs <- vector(mode="list", length=length(fns))
m <- 0
for (i in 1:length(fns)){
    m <- m+1
    print(fns[i])
    TJ2level_CAs[[m]] <- readCA(paste0(fileloc,"/", fns[i]), flexible.symbols = 2)
}
names(TJ2level_CAs) <- TJcat[fns,]$nameInTJ2level_CAs

## check that it works
tjCA(3,21,2)
