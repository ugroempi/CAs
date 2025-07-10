#' Function to make the maximum possible number of rows constant
#'
#' applies maximum clique search for making maximally many rows constant, which
#' may be useful when using an array in a product or power based construction; the constant rows are moved to the front or removed (if requested)
#'
#' @rdname maxconstant
#'
#' @aliases maxconstant
#'
#' @usage maxconstant(D, verbose=0, remove=FALSE, one_is_enough=FALSE, dupcheck=FALSE, ...)
#'
#' @param D a uniform CA without duplicated rows, unless dupcheck is set to TRUE
#' @param verbose 0 for neither comments nor attributes added to output objects,\cr
#' 1 for printing messages with the row numbers used for the constant rows,\cr
#' 2 for not printing anything but adding the list of possible row sets (that could have been made constant) as an attribute to the returned object,\cr
#' 12 for combining both
#' @param remove logical, should constant rows be removed? If TRUE, the function returns a matrix with fewer rows than \code{D}.
#' @param one_is_enough logical; should the algorithm limit itself to making one row constant?
#' @param dupcheck logical, should duplicate rows be checked and removed? If FALSE, the function assumes that there are no duplicate rows.
#' @param ... currently not used
#'
#' @returns an equivalent version of \code{D} that has as many constant rows as possible,
#' located at the start of the array, or removed in case of \code{remove=TRUE}. \cr
#' For \code{verbose=2}, the attribute \code{constant_rows} is attached to the returned object;
#' it is a list with elements \code{design_name} and \code{row_set}, where the latter
#' contains the row numbers in the original design of the rows that were moved to the front.
#'
#' @section Details:
#' The function uses a method presented by Torres-Jimenez at NIST that identifies the maximum
#' number of constant rows via maximum clique size in a graph, whose vertices are the rows
#' of \code{D} and where there are edges between rows that do not coincide in any column.
#' For a clique, all rows can be made constant by the equivalent transform of
#' swapping levels within columns. Of course, the maximum clique size is at most
#' the number of levels, because any further row will have to coincide with an existing one.
#'
#' If the array already has the maximum possible number of constant rows,
#' these are moved to the first runs (or removed), and a search is avoided.\cr
#' If there are as many constant rows in the array as is maximally possible as asserted in a clique search,
#' these are again moved to the first runs (or removed).\cr
#' If there are constant rows that are all together in one of the largest cliques,
#' the first such clique is used.\cr
#' If there are no constant rows and only one constant row is possible, the first row is made constant.\cr
#' For all other cases, the last largest clique is used.\cr
#' The constant rows have ascending levels, starting with the lowest level.
#'
#' @section Warning:
#' For arrays with many rows, the search for the largest cliques can take a long time.
#'
#' @references Torres-Jimenez, Acevedo-Juarez and Avila-George (2021) for the mention of the largest clique idea
#'
#'
#' @examples
#' set.seed(1212) ## for array randomization
#' D <- lhs::createBush(4,5)
#' table(nvals <- lengths(lapply(1:64, function(obj) unique(D[obj,]))))
#' ## one row is already constant:
#' which(nvals==1)
#' ## 4 rows can be made constant (that is the maximum for 4 level CAs)
#' Dmoreconst <- maxconstant(D, verbose=1)
#' head(Dmoreconst)
#' ## previous constant row 22 is now the 2nd row,
#' ## the other clique rows were made constant and moved to positions 1, 3 and 4
#'
#' ## only a single row can be made constant (that is always possible)
#' D <- cyc(19,2) ## k=19, q=19, type 1, t=3
#' table(nvals <- lengths(lapply(1:19, function(obj) unique(D[obj,]))))
#' ## no row is already constant (all have 2 distinct values)
#' ## only a single constant row is achievable
#' head(maxconstant(D, verbose=1))
#'
#' ## two rows can be made constant
#' CA.22.3.12.2 <- cyc(11,2, type="3a") ## N=22, k=12, q=11, type 3a, t=3
#' table(nvals <- lengths(lapply(1:12, function(obj) unique(CA.22.3.12.2[obj,]))))
#' ## no row is already constant (all have 2 distinct values)
#' ## two constant rows are achievable
#' head(Dc <- maxconstant(CA.22.3.12.2, verbose=12))
#' ## which rows of design CA.22.3.12.2 were used ?
#' attributes(Dc)
#'
#' ## example of Colbourn and Torres-Jimenez (2013)
#' ## Figure 2 of the paper
#' hilf <- strsplit(c("2222001", "0022222", "1121201", "0120110",
#'                         "2210202", "2102122", "1200020", "0211121",
#'                         "1012100", "2001210", "0000001", "1111012",
#'                         "2222222", "2222011"), "")
#' plan <- do.call(rbind, lapply(hilf, as.numeric))
#' table(nvals <- lengths(lapply(1:14, function(obj) unique(plan[obj,]))))
#' ## one row is already constant, two rows can be made constant
#' maxconstant(plan)
#'
permvals <- function(D, c, from, to, check=TRUE, ...){
  if (check){
    k <- ncol(D)
    stopifnot(c%%1==0)
    stopifnot(length(from)==length(to))
  }
  if (all(from==to)) return(D) ## no action needed
  hilf <- D[,c]
  for (i in 1:length(from)) D[hilf==from[i],c] <- to[i]
  D
}

#' @importFrom igraph make_empty_graph
#' @importFrom igraph add_edges
#' @importFrom igraph largest_cliques
#' @export
maxconstant <- function(D, verbose=0, remove=FALSE, one_is_enough=FALSE, dupcheck=FALSE, ...){
  Dnam <- deparse(substitute(D))
  stopifnot(verbose %in% c(0,1,2,12))
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  if (dupcheck){
    D <- D[!dupcheck(D),, drop=FALSE]
  }
  N <- nrow(D); k <- ncol(D)
  lev <- levels.no.NA(D)
  if (length(table(lev))>1)
    stop("function maxconstant works for symmetric CAs only")
  v <- nlev <- lev[1]
  consts <- lengths(lapply(1:N, function(obj) unique(D[obj,])))==1
  ## TRUE for constant rows
  nconstant <- sum(consts)

  if (nconstant > 0) posconsts <- which(consts==1)
  if (nconstant > 1) posconsts <- posconsts[sort(D[posconsts,1], index.return=TRUE)$ix]

  ## move the existing maximum possible number of constant rows to the front
  if (nconstant == nlev){
    if (remove)
      D <- D[setdiff(1:N, posconsts),]
    else
      D <- D[c(posconsts, setdiff(1:N, posconsts)),]
    if (verbose %in% c(2,12)) attr(D, "constant_rows") <-
        list(design_name=Dnam, row_set_list=posconsts)
    return(D)
  }

  if (one_is_enough && nconstant==0){
    reflev <- min(D, na.rm=TRUE)
    ## loop over values from second largest to largest
    for (s in (reflev+1):(reflev+v-1)){
      ## treat columns with the value s
      ## setting them to the lowest level
      cs <- which(D[1,]==s)
      if (length(cs) > 0)
        D <- swapvals(D, cs, s, reflev)
    }
    ## now, D has one constant row
    nconstant <- 1
    posconsts <- 1
    ## now handled like the constant row was there already
  }
  if (one_is_enough){
    ## work is finished
    ## more than 1 possible, but only, if already in ingoing D
    nmaxconst <- nconstant
    maxclique <- posconsts
  }else{
    ## further activity is needed
    ## (not yet maximum conceivable number of constant rows,
    ##     and one row was already made constant for one_is_enough)

    ## this should not be needed and might be incomplete, when used in productPCA
    levs <- setdiff(sort(unique(D[,1])), NA) ## removes NA in case of flexible runs
    ## so far, relies on NA values not being in the first few rows

    ## graph with rows as vertices and edges between everywhere distinct rows
    G <- igraph::make_empty_graph(n=N, directed=FALSE)
    paare <- nchoosek(N, 2)
    ## rows of D are vertices
    ## adjacent, if they differ in all columns
    edges <- numeric(0)
    for (i in 1:ncol(paare))
      if (min(abs(D[paare[1,i],]-D[paare[2,i],]), na.rm = TRUE) > 0)
        edges <- c(edges, paare[,i])
    G <- igraph::add_edges(G, edges=edges)

    ## determine all largest cliques
    maxcliques <- igraph::largest_cliques(G)  ## all have same length
    nmaxconst <- length(maxcliques[[1]])

    ## determine the largest clique that will be used
    if (nconstant > 0){
      ## if possible, select a clique that holds all existing constant rows
      candcliques <- which(lengths(lapply(maxcliques, function(obj)
        setdiff(obj, posconsts))) == nmaxconst - nconstant)
      if (length(candcliques)>0) maxclique <- as.numeric(maxcliques[[candcliques[1]]]) else
        maxclique <- as.numeric(rev(maxcliques)[[1]])
    }else{
      ## if only one constant row possible and none is currently constant,
      ## make first row constant
      ## otherwise use the last largest clique
      if (nmaxconst==1) maxclique <- 1 else
        maxclique <- as.numeric(rev(maxcliques)[[1]])  ## rev, because most natural choice usually at the end
    }
    maxclique <- sort(maxclique)
  } ## identification of maxclique finished

  ## general verbosity message
  if (verbose %in% c(1,12))
    message(paste("used largest clique: ", paste0(maxclique, collapse=c(", "))))

  ## no additional constant rows found
  if (nmaxconst==nconstant){
    ## use posconsts, which was resorted above to have rows in ascending order
    ## additional verbosity message
    if (verbose %in% c(1,12)){
      ## additional verbosity message
      message("No additional constant rows found, existing constant rows moved")
      message(paste0("  to rows 1:", nconstant))
    }

    ## output object creation without additional constant rows
    ## but also not v constant rows
    if (remove)
      D <- D[setdiff(1:N, posconsts),]
    else{
      ## place constant rows in front
      D <- D[c(posconsts, setdiff(1:N, posconsts)),]
    }
    ## verbosity attributes for output object
    if (verbose %in% c(2,12)) attr(D, "constant_rows") <-
        list(design_name=Dnam, row_set_list=posconsts)
    return(D)
  }else{
    ## now use maxclique, which potentially contains old constant rows
    ## additional constant row(s) found
    ## move all constant rows to start
    if (remove)
      D <- D[setdiff(1:N, maxclique),]
    else{
      ## rearrange rows
      D <- D[c(maxclique, setdiff(1:N, maxclique)),]
      ## bring first few levels to ascending order
      ## levs is always complete, i.e. nlev==v
      if (nmaxconst == nlev){
        ## can use permvals
        tolevs <- levs[1:nlev]
        for (j in 1:k){
          ## bring all maxclique rows to increasing levels
          D <- permvals(D, j, D[1 : nlev, j], tolevs, check=FALSE)
        }
      }else{
        ## must use swapvals
        for (l in 1:nmaxconst){
          for (j in 1:nlev){
            whichnow <- which(D[l,] == levs[j])
            if (j==l || length(whichnow)==0) next
            D <- swapvals(D, whichnow, levs[j], levs[l]) ## bring first rows to lowest levels
          }
        }
      }}}
  if (verbose %in% c(2,12))
    attr(D, "constant_rows") <-
    list(design_name = Dnam, row_set_list = maxclique)
  D
}

