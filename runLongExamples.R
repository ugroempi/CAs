## run long examples
## \donttest only
devtools::check_built(path="../CAs_0.24-3.tar.gz",
     check_dir = "C:/rtests/ExamInterim", remote=TRUE,
      env_vars = c(NOT_CRAN = "false"))

## \donttest and conditional code
options(run_heavy_examples=TRUE, warn=2)
condexnames <- c("MCA2", "CAEX", "SeqCA_Levenshtein", "optimize_SeqCA")

tryCatch({
for (nam in condexnames){
   message(paste("longrunning examplesIf tests for", nam))
   out <- capture.output({
   print(Sys.time())
   example(nam, character.only=TRUE, run.donttest=TRUE)
   print(Sys.time())
   })
   cat(paste(out, collapse="\n"), file=paste0("output", nam, ".txt"))
}
})
message("Example runs have finished")
options(run_heavy_examples=FALSE, warn=0)
