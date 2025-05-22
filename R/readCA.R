#' Function for reading a covering array from a text file
#'
#' reads a text file into a matrix of class ca (for covering array)
#'
#' @rdname readCA
#'
#' @aliases readCA
#'
#' @usage readCA(path, flexible.symbols=c("*","-","."), comment.symbol="C",
#' ninstruct=1, ignore.chars=NULL, nosep=FALSE, origin=NULL, ...)
#'
#' @param path a character string that specifies the path, including the file name; see Section Details for requirements the file must satisfy.
#' @param flexible.symbols characters to be treated as flexible values (will be set to NA)
#' @param comment.symbol character that starts a comment line
#' @param ninstruct integer, number of lines with instructions (1 or 0)
#' @param ignore.chars \code{NULL}, or characters to be removed from data lines (e.g., \[ and \])
#' @param nosep logical; set to TRUE, if the data lines do not contain separators and each character is a data value
#' @param origin character string to be attached to the output array as origin-attribute;
#' if \code{NULL}, the path and filename information is used
#' @param ... further arguments for function \code{\link{read.table}},
#' which is used separately for each line (after reading them with \code{\link{readLines}})
#'
#' @returns The output object is a matrix of class \code{ca}.
#'
#' @section Details:
#' The file \code{fn} may start with comment lines that begin with a comment symbol, and an instruction line.\cr
#' The first line after removing the comments lines contains instructions about the file content,
#' in the form N k v1^k1 v2^k2 ... It is permissible for the instructions line
#' to have a single string without blanks or \code{^} (e.g., for the strength)
#' after the exponential notation for the columns; this will be ignored in the
#' reading process. \code{ninstruct=0} indicates that there is no instruction line.\cr
#' All columns of the array must have the same starting value, i.e., start all with 0 or all with 1.
#'

#' @importFrom utils read.table
#' @export
readCA <- function(path, flexible.symbols=c("*","-","."), comment.symbol="C",
                   ninstruct=1, ignore.chars=NULL, nosep=FALSE, origin=NULL, ...){
  zeilen <- readLines(con=path)
  zeilen <- zeilen[!nchar(zeilen)==0]
  zeilen <- zeilen[!substr(zeilen,1,1)==comment.symbol]
  stopifnot(ninstruct %in% c(0,1))
  v <- NULL
  if (ninstruct==1){
    ## process instruction information
    instruct <- strsplit(zeilen[1], " ", fixed=TRUE)
    N <- as.numeric(instruct[[1]][[1]])
    k <- as.numeric(instruct[[1]][[2]])
    hilf <- lapply(instruct[[1]][-(1:2)], function(obj)
      as.numeric(unlist(strsplit(obj, "^", fixed=TRUE))))
    hilf <- hilf[lengths(hilf)==2]
    if (length(hilf)==1) uniform <- TRUE else uniform <- FALSE
    if (uniform){
      v <- hilf[[1]][1]
      if (!hilf[[1]][2]==k) stop("contradictory information on k")
    }else{
      ## not uniform
      ## ks from exponential notation
      v <- sapply(hilf, function(obj) obj[1])
      ks <- sapply(hilf, function(obj) obj[2])
      if (!sum(ks)==k) stop("individual ks do not sum to the total k")
    }
    # v <- as.numeric(strsplit(instruct[[1]][[3]],"^", fixed=TRUE)[[1]][[1]])
    zeilen <- zeilen[-1]
    stopifnot(length(zeilen)==N)
  }

  ## zeilen is a character vector
  if (!is.null(ignore.chars)){
    for (ch in ignore.chars){
      zeilen <- gsub(ch, "", zeilen, fixed=TRUE)
    }
  }
  if (nosep){
    zeilen <- funmakefromstrings(zeilen)
    colnames(zeilen) <- paste0("V", 1:ncol(zeilen))
  }
  else{
    zeilen <- t(sapply(zeilen,
                     function(obj) as.matrix(read.table(text=obj,
                                    na.strings=flexible.symbols, ...))))
    rownames(zeilen) <- NULL
  }
  class(zeilen) <- c("ca", class(zeilen))
  if (is.null(origin)) attr(zeilen, "origin") <- path else
    attr(zeilen, "origin") <- origin
  if (is.null(v)) {
    v <- levels.no.NA(zeilen)
    uniform=FALSE
    if (length(table(v))==1){
      uniform <- TRUE
      v <- v[1]
    }
  }
  if (uniform){
     mini <- min(zeilen, na.rm=TRUE);
     maxi <- max(zeilen, na.rm=TRUE)
     stopifnot(maxi-mini==v-1)
  }else{
    # print(ks)
    ## this does not work yet
    # v <- unlist(mapply(rep, v, ks))
      mini <- apply(zeilen, 2, function(obj) min(obj, na.rm=TRUE))
      maxi <- apply(zeilen, 2, function(obj) max(obj, na.rm=TRUE))
     # if (length(unique(mini))>1) stop("all columns must have the same minimum value")
    #  if (!all(maxi - mini == v - 1)) stop("there is a problem with the coding of the columns")
  }
  zeilen
}

