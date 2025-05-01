#' Function for reading a covering array from a text file
#'
#' reads a text file into a matrix of class ca (for covering array)
#'
#' @rdname readCA
#'
#' @aliases readCA
#'
#' @usage readCA(path, fn, flexible.symbols=c("*","-","."), comment.symbol="C",
#' origin=NULL, ...)
#'
#' @param path a character string that specifies the path
#' @param fn a character string that specifies the file name; see Section Details for requirements the file must satisfy.
#' @param flexible.symbols characters to be treated as flexible values (will be set to NA)
#' @param comment.symbol character that starts a comment line
#' @param origin character string to be attached to the output array as origin-attribute;
#' if \code{NULL}, the path and filename information is used
#' @param ... further arguments for function \code{\link{read.table}},
#' which is used separately for each line (after reading them with \code{\link{readLines}})
#'
#' @returns The output object is a matrix of class \code{ca}.
#'
#' @section Details:
#' The file \code{fn} must start with comment lines that begin with a comment symbol.\cr
#' The first line after the comments lines contains instructions about the file content,
#' in the form N k v1^k1 v2^k2 ... It is permissible for the instructions line
#' to have a single string without blanks or \code{^} (e.g., for the strength)
#' after the exponential notation for the columns; this will be ignored in the
#' reading process.\cr
#' All columns of the file must have the same starting value, i.e., start all with 0 or all with 1.
#'

#' @importFrom utils read.table
#' @export
readCA <- function(path, fn, flexible.symbols=c("*","-","."), comment.symbol="C", origin=NULL, ...){
  zeilen <- readLines(con=paste0(path, "\\", fn))
  zeilen <- zeilen[!nchar(zeilen)==0]
  zeilen <- zeilen[!substr(zeilen,1,1)==comment.symbol]
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
    v <- sapply(hilf, function(obj) obj[1])
    ks <- sapply(hilf, function(obj) obj[2])
    if (!sum(ks)==k) stop("individual ks do not sum to the total k")
  }
  # v <- as.numeric(strsplit(instruct[[1]][[3]],"^", fixed=TRUE)[[1]][[1]])
  zeilen <- zeilen[-1]
  stopifnot(length(zeilen)==N)

  ## zeilen is a character vector
  zeilen <- t(sapply(zeilen,
                     function(obj) as.matrix(read.table(text=obj,
                                    na.strings=flexible.symbols, ...))))
  rownames(zeilen) <- NULL
  class(zeilen) <- c("ca", class(zeilen))
  pfad <- paste0(path, "\\", fn)
  if (is.null(origin)) attr(zeilen, "origin") <- pfad else
    attr(zeilen, "origin") <- origin
  if (uniform){
     mini <- min(zeilen, na.rm=TRUE);
     maxi <- max(zeilen, na.rm=TRUE)
     stopifnot(maxi-mini==v-1)
  }else
  {
    v <- unlist(mapply(rep, v, ks))
      mini <- apply(zeilen, 2, function(obj) min(obj, na.rm=TRUE))
      maxi <- apply(zeilen, 2, function(obj) max(obj, na.rm=TRUE))
      if (length(unique(mini))>1) stop("all columns must have the same minimum value")
      if (!all(maxi - mini == v - 1)) stop("there is a problem with the coding of the columns")
  }
  zeilen
}

