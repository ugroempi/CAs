#' Function for reading a covering array from a text file
#'
#' reads a text file into a matrix of class ca (for covering array)
#'
#' @rdname readCA
#'
#' @aliases readCA
#'
#' @usage readCA(path, flexible.symbols=c("*","-","."), comment.symbol="C",
#' ninstruct=1, skiplines=0, ignore.chars=NULL, nosep=FALSE, origin=NULL, ...)
#'
#' @param path a character string that specifies the path, including the file name; see Section Details for requirements the file must satisfy.
#' @param flexible.symbols characters to be treated as flexible values (will be set to NA)
#' @param comment.symbol character that starts a comment line
#' @param ninstruct integer, number of lines with instructions (1 or 0)
#' @param skiplines integer, number of lines with no readable instructions and no defined comment character
#'        that should nevertheless be skipped; those are assumed to be above the instructions line
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
#' The file may contain blank lines and comment lines whose first non-whitespace
#' character is the comment symbol. These lines are ignored.\cr
#' The first line after removing the comments lines contains instructions about the file content,
#' in the form N k v1^k1 v2^k2 ..., with fields separated by spaces or tabs.
#' It is permissible for the instructions line
#' to have a single string without blanks or \code{^} (e.g., for the strength)
#' after the exponential notation for the columns; this will be ignored in the
#' reading process. \code{ninstruct=0} indicates that there is no instruction line.\cr
#' Each data line must contain the same number of entries. When an instruction
#' line is present, the data must have exactly N rows and k columns. A single
#' row or column is retained as a matrix. These checks validate the file shape,
#' not the covering strength of the array.\cr
#' All columns of the array must have the same starting value, i.e., start all with 0 or all with 1.
#'

#' @importFrom utils read.table
#' @export
readCA <- function(path, flexible.symbols=c("*","-","."), comment.symbol="C",
                   ninstruct=1, skiplines=0, ignore.chars=NULL, nosep=FALSE, origin=NULL, ...){
  zeilen <- readLines(con=path)
  zeilen <- zeilen[nzchar(trimws(zeilen))]
  zeilen <- zeilen[!substr(trimws(zeilen),1,1)==comment.symbol]
  if (skiplines > 0) zeilen <- zeilen[-(1:skiplines)]
  stopifnot(ninstruct %in% c(0,1))
  if (length(zeilen)==0) stop("no array data found")
  v <- NULL
  if (ninstruct==1){
    ## process instruction information
    instruct <- strsplit(trimws(zeilen[1]), "[[:space:]]+")
    dims <- suppressWarnings(as.numeric(instruct[[1]][1:2]))
    if (anyNA(dims) || any(!is.finite(dims)) ||
        any(dims < 1 | dims %% 1 != 0))
      stop("instruction line must start with positive integer N and k")
    N <- dims[1]
    k <- dims[2]
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
    if (length(zeilen)!=N)
      stop("number of data rows does not match N in the instruction line")
  }

  ## zeilen is a character vector
  if (!is.null(ignore.chars)){
    for (ch in ignore.chars){
      zeilen <- gsub(ch, "", zeilen, fixed=TRUE)
    }
  }
  if (nosep){
    if (length(unique(nchar(zeilen)))!=1)
      stop("data rows must all have the same number of entries")
    zeilen <- funmakefromstrings(zeilen)
    colnames(zeilen) <- paste0("V", 1:ncol(zeilen))
  }
  else{
    rows <- lapply(zeilen, function(obj) as.matrix(read.table(text=obj,
                                    na.strings=flexible.symbols, ...)))
    widths <- vapply(rows, ncol, integer(1))
    if (any(vapply(rows, nrow, integer(1))!=1) ||
        any(widths==0) || length(unique(widths))!=1)
      stop("data rows must all have the same number of entries")
    zeilen <- do.call(rbind, rows)
    dimnames(zeilen) <- list(NULL, NULL)
  }
  if (ninstruct==1 && ncol(zeilen)!=k)
    stop("number of data columns does not match k in the instruction line")
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

