library(CAs)

# Test file parsing and matrix shape, independently of covering-strength checks.
check_error <- function(expr, pattern) {
  result <- tryCatch(force(expr), error=identity)
  stopifnot(inherits(result, "error"),
            grepl(pattern, conditionMessage(result), fixed=TRUE))
}

# readCA accepts whitespace in headers, comments and blank lines
local({
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(c("  C source comment", "\t ", "  4\t 3   2^3  2 ",
               "0 0 0", "0 1 1", "  ", "1 0 1", "1 1 0"), path)
  D <- readCA(path)
  stopifnot(inherits(D, "ca"))
  stopifnot(identical(dim(D), c(4L, 3L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0))))
  stopifnot(identical(attr(D, "origin"), path))
})

# readCA retains single rows, single columns and flexible entries
local({
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(c("1 3 2^3", "0 1 *"), path)
  D <- readCA(path, origin="one row")
  stopifnot(identical(dim(D), c(1L, 3L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 1, NA))))
  stopifnot(identical(attr(D, "origin"), "one row"))

  writeLines(c("3 1 2^1", "0", "1", "*"), path)
  D <- readCA(path)
  stopifnot(identical(dim(D), c(3L, 1L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 1, NA))))

  writeLines(c("0,1,*", "1,0,1"), path)
  D <- readCA(path, ninstruct=0, sep=",", header=FALSE)
  stopifnot(identical(dim(D), c(2L, 3L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 1, 1, 0, NA, 1))))
})

# readCA rejects data that disagree with the declared shape
local({
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(c("4 3 2^3", "0 0", "0 1", "1 0", "1 1"), path)
  check_error(readCA(path), "number of data columns does not match k")
  writeLines(c("3 2 2^2", "0 0", "1 1"), path)
  check_error(readCA(path), "number of data rows does not match N")
  writeLines(c("4 x 2^3", "0 0 0"), path)
  check_error(readCA(path), "positive integer N and k")
  writeLines(c("C comments only", "  "), path)
  check_error(readCA(path), "no array data found")
})

# readCA rejects ragged rows with and without separators
local({
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(c("0 1", "1"), path)
  check_error(readCA(path, ninstruct=0), "same number of entries")
  writeLines(c("01", "1"), path)
  check_error(readCA(path, ninstruct=0, nosep=TRUE), "same number of entries")
  writeLines(c("2 3 2^3", "01", "10"), path)
  check_error(readCA(path, nosep=TRUE), "number of data columns does not match k")
})

# readCA preserves bracketed, compact and mixed-level formats
local({
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(c("[0,1]", "[1,0]"), path)
  D <- readCA(path, ninstruct=0, ignore.chars=c("[", "]"),
              sep=",", header=FALSE)
  stopifnot(identical(dim(D), c(2L, 2L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 1, 1, 0))))

  writeLines(c("2 2 2^2", "01", "10"), path)
  D <- readCA(path, nosep=TRUE)
  stopifnot(identical(dim(D), c(2L, 2L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 1, 1, 0))))
  stopifnot(identical(colnames(D), c("V1", "V2")))

  writeLines(c("ignored heading", "3 2 2^1 3^1", "0 0", "1 1", "0 2"), path)
  D <- readCA(path, skiplines=1)
  stopifnot(identical(dim(D), c(3L, 2L)))
  stopifnot(isTRUE(all.equal(as.vector(D), c(0, 1, 0, 0, 1, 2))))
})
