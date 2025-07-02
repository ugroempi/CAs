nist_folder <- "https://math.nist.gov/coveringarrays/ipof/cas"
## it is not permitted to download directory listings from there
## https://math.nist.gov/coveringarrays/ipof/ipof-results.html
## https://math.nist.gov/coveringarrays/ipof/tables/table.2.2.html etc.
## example file - these can be downloaded separately
##    can be automated, but requires a lot of manual intervention
## https://math.nist.gov/coveringarrays/ipof/cas/t=2/v=2/ca.2.2^3.txt.zip

## NISTcatpath <- path where the downloaded zip files reside

zipfiles <- list.files(path=NISTcatpath)
## remove unsuitable files from listing
zipfiles <- setdiff(zipfiles, c("download_NIST_CAs.R", "NISTcat.RData", "colbourn.R"))

ts <- sapply(zipfiles, function(obj) as.numeric(substr(obj,4,4)))
vs <- sapply(zipfiles, function(obj) as.numeric(substr(obj,6,6)))
ks <- sapply(zipfiles, function(obj) as.numeric(substr(obj,8, regexpr(".txt.zip", obj, fixed=TRUE)-1)))

## run sizes must be determined from file content
## runs several minutes
system.time(
sizes <- sapply(1:length(ts), function(obj) readr::read_table(file=paste0(NISTcatpath, "\\", "ca.", ts[obj],".",vs[obj],"^",ks[obj],".txt.zip"), n_max=1, col_names=FALSE, col_types="i"))
)
#       User      System     elapsed 
#      97.75       38.75      289.75 

NISTcat <- cbind(t=ts, v=vs, k=ks, N=unlist(sizes))