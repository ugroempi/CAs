#' function to make the maximumm possible number of rows constant
#'
#' applies maximum clique search for making maximally many rows constant, which
#' may be useful when using an array in a product or power based construction; the constant rows are moved to the front or removed (if requested)
#'
#' @rdname maxconstant
#'
#' @aliases maxconstant
#'
#' @usage maxconstant(D, verbose=0, remove=FALSE, dupcheck=FALSE, ...)
#'
#' @param D a uniform CA without duplicated rows, unless dupcheck is set to TRUE
#' @param verbose 0 for no comments, 1 for printing messages with the row numbers used for the constant rows,
#' 2 for additionally adding the list of possible row sets (that could have been made constant) as an attribute to the returned object
#' @param remove logical, should constant rows be removed? If TRUE, the function returns a matrix with fewer rows than \code{D}.
#' @param dupcheck logical, should duplicate rows be checked and removed? If FALSE, the function assumes that there are no duplicate rows.
#' @param ... currently not used
#'
#' @returns an equivalent version of \code{D} that has as many constant rows as possible,
#' located at the start of the array, or removed in case of \code{remove=TRUE}. For \code{verbose=2}, a list of the design name and the list of possible rows
#' to be made constant is attached as the attribute \code{potential_constant_rows} to the returned object
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
#' For all other cases, the first largest clique is used.\cr
#' Making rows constant changes all levels to those of the first column.
#'
#' @section Warning:
#' For arrays with many rows, the search for the largest cliques can take a long time.
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
#' head(Dc <- maxconstant(CA.22.3.12.2, verbose=2))
#' ## which rows of design CA.22.3.12.2 could have been used (including the ones that were used)?
#' attributes(Dc)
#'
#' D <- KSK(k=12)
#' maxconstant(D, verbose=1)
#'
permvals <- function(D,c,from,to,check=TRUE,...){
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
maxconstant <- function(D, verbose=0, remove=FALSE, dupcheck=FALSE, ...){
   Dnam <- deparse(substitute(D))
   stopifnot(verbose %in% c(0,1,2))
   stopifnot(is.matrix(D))
   stopifnot(is.numeric(D))
   if (dupcheck){
     D <- D[!dupcheck(D),, drop=FALSE]
   }
   N <- nrow(D); k <- ncol(D)
   lev <- levels.no(D)
   if (length(table(lev))>1)
     stop("function maxconstant works for symmetric CAs only")
   nlev <- lev[1]
   consts <- lengths(lapply(1:N, function(obj) unique(D[obj,])))==1
   ## TRUE for constant rows
   nconstant <- sum(consts)
   if (nconstant > 0) posconsts <- which(consts==1)
   ## move the existing maximum possible number of constant rows to the front
   if (nconstant == nlev){
     if (remove)
       D <- D[setdiff(1:N, posconsts),]
     else
       D <- D[c(posconsts, setdiff(1:N, posconsts)),]
   }else{
     ## not yet maximum conceivable number of constant rows
   levs <- sort(unique(D[,1]))
   G <- igraph::make_empty_graph(n=N, directed=FALSE)
   paare <- nchoosek(N, 2)
   ## rows of D are vertices
   ## adjacent, if they differ in all columns
   edges <- numeric(0)
   for (i in 1:ncol(paare))
     if (min(abs(D[paare[1,i],]-D[paare[2,i],])) > 0)
         edges <- c(edges, paare[,i])
   G <- igraph::add_edges(G, edges=edges)
   ## determine all largest cliques
   maxcliques <- igraph::largest_cliques(G)
   nmaxconst <- length(maxcliques[[1]])
   if (nconstant>0){
     ## determine whether any clique holds all existing constant rows
     candcliques <- which(lengths(lapply(maxcliques, function(obj)
       setdiff(obj, posconsts)))==nmaxconst-nconstant)
     if (length(candcliques)>0) maxclique <- unclass(maxcliques[[candcliques[1]]]) else
       maxclique <- unclass(maxcliques[[1]])
   }else {
     ## if only one constant row possible and none is currently constant,
     ## make first row constant
     ## otherwise use the first largest clique
     if (nmaxconst==1) maxclique <- 1 else
     maxclique <- as.numeric(maxcliques[[1]])
   }
   if (verbose>0)
     message(paste("used largest clique: ", paste0(maxclique, collapse=c(", "))))
   if (length(maxclique)==nconstant){
     ## use the existing constant rows
     if (verbose>0) {
       message("No additional constant rows found, existing constant rows moved")
       message(paste0("  to rows 1:", nconstant))
     }
     if (remove)
       D <- D[setdiff(1:N, posconsts),]
     else
       D <- D[c(posconsts, setdiff(1:N, posconsts)),]
   }else{## additional constant row(s) found
     if (length(maxclique)>1){
     tolevs <- D[maxclique, 1] ## bring all other columns to this combination
     for (j in 2:k){
       ## bring all maxclique rows to constant value of first column
       D <- permvals(D, j, D[maxclique,j], tolevs, check=FALSE)
     }
     ## move constant rows to start
     if (remove)
       D <- D[setdiff(1:N, maxclique),]
     else
       D <- D[c(maxclique,setdiff(1:N, maxclique)),]
     }else{ ## length of maxclique is = 1
       for (j in 2:k){
         D <- swapvals(D,j,D[1,1],D[1,j])
       }}
   }}
 if (verbose==2) {
   if (nconstant==nlev) attr(D, "possible_constant_rows") <- list(design_name=Dnam, row_set_list=maxclique) else
   attr(D, "possible_constant_rows") <- list(design_name=Dnam, row_set_list=maxcliques)
 }
 D
}
