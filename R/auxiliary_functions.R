#' auxiliary functions
#'
#' non-visible functions to be used in other functions,
#' and currently also previous drafts that may be deleted
#'
#' @rdname auxiliary
#'
#' @usage NULL
#'
#' @details
#'
#' matcheck: checks a matrix for suitability
#' gf_minus: Galois field subtraction
#' funmakefromstrings: read example designs from copy-pasted strings
#' iacheck: check a two-column matrix for correct interaction structure
#' subia: checks whether ia1 is subset of ia2
#' levels.no from DoE.base
#' fasttab: table replacement without using factors
#' LS_to_D3cols: function for making a latin square into a three column design
#' D3cols_to_LS: function for bringing a three column design to latin square format
#' metrics: a function for calculation metrics, likely obsolete
#' proportions_slower: earlier version of coverage, likely obsolete
#' ord from DoE.base
#' rho: extracts run numbers of D that hold the interaction ia
#' allints: function that extracts all interactions for (up to) t factors from a run
#' failcands: fail candidates from a test suite with a single failed run
#'   failcands2: (failed) attempt to speed up failcands
#' rhos_for_ints
#' is.CA: checks for covering array properties,
#'   can return distribution of covering frequencies or actual row sets
#' is.LA: checks for locating array properties
#' is.DA: checks for detecting array properties (not ready)
#'
#' @importFrom utils tail

## Function definition
## from https://stackoverflow.com/questions/17683370/how-to-produce-combinations-iteratively-in-r
gen.next.cbn <- function(cbn, n){
  ## Generates the combination that follows the one provided as input
  cbn.bin      <- rep(0, n)
  cbn.bin[cbn] <- 1
  if (tail(cbn.bin, 1) == 0){
    ind <- tail(which(cbn.bin == 1), 1)
    cbn.bin[c(ind, ind+1)] <- c(0, 1)
  }else{
    ind <- 1 + tail(which(diff(cbn.bin) == -1), 1)
    nb  <- sum(cbn.bin[-c(1:ind)] == 1)
    cbn.bin[c(ind-1, (n-nb+1):n)] <- 0
    cbn.bin[ind:(ind+nb)]         <- 1
  }
  cbn <- which(cbn.bin == 1)
}

### Iteration example
### Example parameters
#n   <- 6
#k   <- 3
#
#system.time({
#for (i in 1:choose(n, k)){
#  if (i == 1){
#    cbn <- 1:k
#  }else{
#    cbn <- gen.next.cbn(cbn, n)
#  }
#  #print(cbn)
#}
#})


matcheck <- function(mat, start0=NULL, uniform=TRUE,
                     PCAcheck=FALSE, flexible=NULL, ...){
  ## PCAcheck checks for the attribute
  ##    (a check for a PCA structure is done with is.PCA, but only picking up on the unchanged array)
  ## flexible=NULL does check for NA values and treats these as flexible
  ##    other values must be declared

  ## flexible (for the don't care value symbol) must be NA or numeric, I suspect that it is v for Torres-Jimenez
  stopifnot(is.matrix(mat))
  stopifnot(is.numeric(mat))
  stopifnot(all(mat%%1==0, na.rm = TRUE))

  ## flexible is taken from an attribute, if present
  if (is.null(flexible)){
    if (!is.null(flexibleAttr <- attr(mat,"flexible")))
      flexible <- flexibleAttr$value
  }
  if (is.null(flexible)){
    if (any(is.na(mat))) flexible <- NA
  }
  if (!is.null(flexible)){
    ## flexible values are fixed for the rest of the check
    stopifnot(is.na(flexible) || is.numeric(flexible))
    stopifnot(length(flexible)==1)
    anyflex <- FALSE
        ## fixed level
           if (is.na(flexible)){
               if (any(is.na(mat))){
                 anyflex <- TRUE
                 mat[is.na(mat)] <- 1
               }
             }else
               if (any(mat==flexible)){
                 anyflex <- TRUE
                 mat[mat==flexible] <- 1
               }
           ## still length(flexible)==1
           if (!anyflex) flexible <- NULL
    ## for DECLARED flexible values:
    ## now there is a valid entry instead of NA (or whatever else),
    ## flexible has been set NULL, if no such values occur
    ## flexible mat values have been replaced with 1
  }
  if (is.null(start0)){
       stopifnot(min(mat) %in% c(0,1))
       if (all(apply(mat,2,function(obj) min(obj))==0))
         start0 <- TRUE else
           if (all(apply(mat,2,min)==1)) start0 <- FALSE
    }else{
    if (start0) stopifnot(min(mat)==0) else
      stopifnot(min(mat)==1)
  }
  ll <- levels.no(mat)
  if (PCAcheck){
    PCAstatus <- attr(mat, "type")
    if (is.null(PCAstatus)) {
      PCAstatus <- is.PCA(mat)
      if (PCAstatus) PCAstatus <- attr(PCAstatus, "PCAstatus")
      ## otherwise FALSE anyway
      }
  }
  if (uniform){
      stopifnot(length(unique(ll))==1)
    aus <- list(v=ll[1], k=ncol(mat), N=nrow(mat))
  } else{
    aus <- list(vs=ll, k=ncol(mat), N=nrow(mat))
  }
  if (PCAcheck) aus <- c(aus, PCAstatus=list(PCAstatus))
  if (!is.null(flexible)) aus <- c(aus, flexvalue=flexible)
  if (is.null(start0)) aus <- c(aus, start0="mixed") else
                       aus <- c(aus, start0=start0)
  aus
}

cycvec <- function(v, q, gf=NULL, primitive=NULL){
  ## t is the requested coverage
  ## v is the common number of levels
  ## q is the prime or prime power on which the construction is based
  if (is.null(gf)) gf <- lhs::create_galois_field(q)
  ## start vector is the vector of "logarithms of 1:q" (mod v) <w.r.t. the
  ##     base chosen as a primitive element (omega) of the group based on q
  ## retrieve primitive
  if (!is.null(primitive)) p <- primitive else p <- primedat[which(primedat[,1]==q),2]
  ## create start vector
  if (!q %% v==1) stop("q mod v = 1 is violated")

  xstart <- rep(0,q)
  cur <- p  ## current primitive
  for (i in 1:(q-1)) {
    xstart[cur+1] <- i  ## exponent of the primitive is 1
    cur <- SOAs:::gf_prod(cur, p, gf) # next exponent
  }
  if (!all(table(xstart)==1)) stop("wrong primitive")
  xstart%%v
}

gf_minus <- function(x,y,gf){
  ## calculates x - y, after reducing both by a mod operation to
  ## gf entries
  q <- gf$q; p <- gf$p
  x <- x%%q; y <- y%%q
  yn <- gf$neg[y+1]
  gf$plus[x+1,yn+1]
}


## function for reading and arranging copy-pasted strings
funmakefromstrings <- function(stringvec){
  do.call(rbind, lapply(strsplit(stringvec, "", fixed=TRUE), as.numeric))
}

#' @exportS3Method NULL
#'
levels.no.NA <- function (x){
  xx <- x
  ff <- FALSE
  if (is.data.frame(xx)) {
    if (any(ff <- sapply(xx, is.factor)))
      nflevs <- sapply(xx[ff], nlevels)
  }
  aus <- apply(xx, 2, function(v) length(setdiff(unique(v),NA)))
  if (any(ff))
    aus[ff] <- nflevs
  aus
}

iacheck <- function(ia, t=NULL, lno=NULL, upto=FALSE, start0=TRUE){
  ## two column matrix with first index a factor,
  ## second index its value
  ## first index must not have duplicates!
  ##    needs check for valid ia format (iacheck)
  ##    needs subset method between interactions (subia)
  ##    needs row selection method of interaction and design (rho)
  if (!is.matrix(ia)) return(FALSE)
  if (is.null(t)) t <- nrow(ia)
  ## wrong dimension
  if (!(nrow(ia)==t || (upto && nrow(ia)<t)) || !ncol(ia)==2) return(FALSE)
  ## non-unique ia[,1]
  if (!length(unique(ia[,1]))==nrow(ia)) return(FALSE)
  ## non-integer levels
  if (!all(ia[,2]%%1==0)) return(FALSE)
  ## non-permissible levels
  if (!all(ia[,2]>=ifelse(start0, 0, 1))) return(FALSE)
  ## invalid levels
  if (!is.null(lno))
  if (!all(lno[ia[,1]]-start0 >= ia[,2])) return(FALSE)
  TRUE
}

#iacheck(ia1, 2)

subia <- function(ia1, ia2){
  ## checks whether ia1 is a subset of ia2
  t1 <- nrow(ia1); t2 <- nrow(ia2)
  stopifnot(iacheck(ia1, t1) && iacheck(ia2, t2))
  if (t1 > t2) return(FALSE)
  if (t1==t2)
    if (all(ia1[ord(ia1),]==ia2[ord(ia2),])) return(TRUE) else return(FALSE)
  ## now nrow(ia1) < nrow(ia2)
  if (length(andere <- setdiff(ia2[,1],ia1[,1])) > t2-t1) return(FALSE)
  ia2compare <- ia2[-which(ia2[,1]%in%andere),]
  ia1 <- ia1[ord(ia1),]
  ia2compare <- ia2[ord(ia2compare),]
  all(ia1==ia2compare)
}

#' @exportS3Method NULL
#'
levels.no <- DoE.base:::levels.no
ord <- DoE.base::ord
nchoosek <- DoE.base:::nchoosek
## check for rows with the interaction in design D
rho <- function(D, ia, start0=TRUE){
  ## D the design
  ## ia the interaction (tx2 matrix)
  ## start0 levels start with 0 (otherwise with 1)
  lno <- levels.no(D)
  stopifnot(iacheck(ia, nrow(ia), lno, start0=start0))
  stopifnot(max(ia[,1])<=ncol(D))
  t <- nrow(ia)
  ## actual values of the columns of interest (as a list of column vectors)
  hilf <- lapply(ia[,1], function(obj) D[,obj])
  ## function to return row numbers
  ##    for which the respective values are taken on
  myfun <- function(ll1, ll2)
    lapply(1:length(ll1), function(obj) which(ll1[[obj]]==ll2[[obj]]))
  ## count the row numbers,
  ## return row numbers for which all t coincidences are met
  hilf <- table(do.call(c, myfun(hilf, as.list(ia[,2]))))
  as.numeric(names(hilf)[hilf==t])
}
#rho(d2f, ia1)
#rho(d2f, ia2)
#rho(d1d, ia2)
#rho(d1e, ia2)

allints <- function(D, r, t, upto=FALSE){
  # D design
  # r run number or vector of run numbers
  # t strength (t-way interactions)
  # upto: if TRUE, up to t-way interactions
  # returns list of interactions
  if (upto && t>1){
    ll <- rep(list(NULL), t)
    for (i in 1:t){
      ll[[i]] <- allints(D,r,i,upto=FALSE)
    }
    ll <- unlist(ll, recursive=FALSE)
  }else{
    hilf <- nchoosek(ncol(D),t)
    if (length(r)==1)
    ll <- lapply(1:ncol(hilf),
                 function(obj) cbind(hilf[,obj], D[r,hilf[,obj]]))
    else
      ll <- unlist(lapply(1:ncol(hilf),
                          function(obj)
                            lapply(r, function(obj2)
                              cbind(hilf[,obj], D[obj2,hilf[,obj]]))),
                          recursive=FALSE)
  }
  ll
}

failcands <- function(D, r, t, upto=FALSE){
  # D a test suite (design)
  # r a run number or a vector of run numbers of failed tests
  # t the degree of interactions
  # upto logical that says whether interactions with up to t or exactly t factors

  ## fail candidates from a single failed run
  ## at least one must be actually a failure cause

  ## fail candidates from a vector of failed runs
  ## at least one must be a failure cause, likely more than one needed,
  ##     because all the runs have failed
  candfail <- allints(D, r, t, upto=upto)
  n <- nrow(D)
  candpass <- unique(unlist(lapply(setdiff(1:n,r),
                                   function(obj) allints(D, obj, t, upto=upto)),
                            recursive=FALSE)
                     )
  setdiff(candfail, candpass)
}

failcands2 <- function(D, r, t, upto=FALSE){
  ## fail candidates from a single failed run
  ## works for upto=FALSE only, at the moment
  ## and is not faster but slower
  ## nevertheless, might be necessary for large numbers of columns
  candfail <- allints(D, r, t, upto=upto)
  facs <- 1:ncol(D)
  andere <- setdiff(1:nrow(D), r)
  for (i in andere){
    weg <- allints(D[i,facs,drop=FALSE], 1, t, upto=upto)
    weg <- lapply(weg, function(obj) {
      obj[,1] <- facs[obj[,1]]
      obj
    })
    candfail <- setdiff(candfail, weg)
    facs <- sort(unique(c(sapply(candfail, function(obj) obj[,1]))))
  }
  candfail
}

rhos_for_ints <- function(D, t, upto=FALSE, start0=TRUE){
  ## calculates the sets of rows for each interaction
  lnos <- levels.no.NA(D)
  m <- ncol(D)
  ints_per_row <- choose(m, t)
  pickcols <- nchoosek(m,t)
  liste <- vector(mode="list", length=ints_per_row)
  names(liste) <- sapply(1:ncol(pickcols), function(obj)
    paste(pickcols[,obj], collapse="-"))
  if (start0) {
    for (i in 1:ncol(pickcols)){
      hilf <- lapply(lnos[pickcols[,i]], function(obj) 0:(obj-1))
      allintsi <- as.matrix(do.call(expand.grid, hilf))
        ## has only the relevant columns
      hilf <- lapply(1:nrow(allintsi), function(obj)
        rho(D, cbind(pickcols[,i], allintsi[obj,]), start0 =TRUE))
      names(hilf) <- apply(allintsi, 1, function(obj) paste(obj, collapse="."))
      liste[[i]] <- hilf
    }
  } else{
    for (i in 1:ncol(pickcols)){
      hilf <- lapply(lnos[pickcols[,i]], function(obj) 1:obj)
      allintsi <- as.matrix(do.call(expand.grid, hilf))
        ## has only the relevant columns
      hilf <- lapply(1:nrow(allintsi), function(obj)
        rho(D, cbind(pickcols[,i], allintsi[obj,]), start0 = FALSE))
      names(hilf) <- apply(allintsi, 1, function(obj) paste(obj, collapse="."))
      liste[[i]] <- hilf
    }
  }
  ## important: all elements of liste are always sorted in increasing order
  return(liste)
}

is.CA <- function(D, t, index=1, start0=TRUE, verbose=0){
  stopifnot(is.logical(start0))
  allrhos <- rhos_for_ints(D, t, upto=TRUE, start0=start0)
  if (verbose==0) return(all(lengths(unlist(allrhos, recursive=FALSE))>=index))
  aus <- all(lengths(unlist(allrhos, recursive=FALSE))>0)
  if (verbose==2) attr(aus, "rowsets") <- allrhos
  if (verbose==1){
    attr(aus, "frequDistr") <- table(lengths(unlist(allrhos, recursive=FALSE)))
    attr(aus, "proportions") <- sapply(allrhos, function(obj) sum(lengths(obj)>0)/length(obj))
  }
  aus
}

is.LA <- function(D, t, d, upto=FALSE, start0=TRUE, dupto=FALSE, verbose=0){
  ## implemented variants: d=1 should be correct
  ## d=2 and t=2 not yet trustworthy
  ## d=2 and t=1 not yet finished
  stopifnot(is.logical(upto))
  stopifnot(is.logical(dupto))
  stopifnot(is.logical(start0))
  if (!is.CA(D,t,start0=start0)) return(FALSE)
  m <- ncol(D)
  ## rho sets for individual t-factor interactions
  allrhos <- rhos_for_ints(D, t, upto=upto, start0=start0)
  ## a list-valued entry for each t-factor combination
  ##     length of the list element is the number of possible instances
  ##     of interactions for the specific set of t factors

  ## dups holds TRUE in the element positions that are identical
  ##    to earlier elements of the list of row number vectors,
  ##    FALSE for first occurrence or unique values
  dups <- duplicated(unlist(allrhos, recursive=FALSE))
  ##################
  ## d=1, dupto irrelevant, upto (for t) handled in allrhos
  if (d==1) {
    if (verbose==0) return(!any(dups))
    else{
      aus <- !any(dups)
      attr(aus, "liste") <- allrhos
      return(aus)
    }
  }
  ##################

  ##################
  if (d==2 && t==1){
    cases <- nchoosek(m, t)  ## entries are column numbers of D
                             ## 1 x m matrix, entries 1:m
    sammel <- array(data=vector(mode="list"), dim=c(m, m))
    rownames(sammel) <- colnames(sammel) <- cases
    anzahlen <- lengths(allrhos) ## coincides with levels.no(D)
    ## lengths of single vectors (for upto=TRUE)

    ## index number of levels for the columns of d (starting with 1)
    ## equal to index number of the 1-way interactions
    vecs <- mapply(function(obj) 1:obj, anzahlen)

    ## combinations of d different interactions
    dset <- nchoosek(length(anzahlen), d)
    ## add cases with same columns
    ## for levels.no>2 ## (otherwise would always be all rows)
    nlevs <- levels.no(D)  ## same as anzahlen
    hilf <- (1:m)          ## for more than one factor with 2 levels
                           ## LA(2,1) is never possible
    if (length(hilf)>0)
      dset <- cbind(dset, rbind(hilf,hilf))
    ## later perhaps handle with gtools::combinations with repeats.allowed = TRUE
    nintcombis <- ncol(dset)
    for (i in 1:nintcombis){
      cols <- dset[,i]  ## vector of length d with d vector combination indices in it
      {
        x <- do.call(expand.grid, vecs[cols])
        ## the selected combinations
        ## has d=2 columns (for the d vector combination indices)
        ## same vector combinations must use different instances
        ## and each instance only once
        if (cols[1]==cols[2]){
          if (dupto)
          x <- x[x[,1]<=x[,2],] ## including = for d=1
          else
          x <- x[x[,1] < x[,2],]
        }
        sammel[min(cols), max(cols)][[1]] <-
          lapply(1:nrow(x), ## every combination that was not removed
                 function(obj){
                   ## obj is the row no of x
                   ## obj2 is the list of rho sets for the tfis of dset[,i],
                   ## j is the row of x (i.e. selected interaction instance)
                   sort(Reduce(union,
                          mapply(function(obj2, j) obj2[[j]],
                                 obj2=allrhos[dset[,i]],
                                 j=x[obj,],
                                 SIMPLIFY = FALSE)))
                 })
        ## elements of sammel hold the list of row sets for unions of two instances
        ##     of combined interactions
        if (verbose>0) names(sammel[min(cols), max(cols)][[1]]) <-
          apply(as.matrix(x), 1, function(obj) paste(obj, collapse=":"))
      }
    }
    hilf <- c(sammel)  ## array structure removed
    if (verbose>0) {
      hilf2 <- expand.grid(rownames(sammel), colnames(sammel))
      names(hilf) <- paste(hilf2[,1], hilf2[,2], sep="-")
    }
    hilf <- hilf[sapply(hilf, function(obj) !is.null(obj))]
    if (dupto){
      ## d=1
      allrhos <- rhos_for_ints(D, t, upto=upto, start0=start0)
      ## a list-valued entry for each t-factor combination
      ##     length of the list element is the number of possible instances
      ##     of interactions for the specific set of t factors
      hilf <- c(unlist(allrhos, recursive=FALSE), hilf)
    }
    ## t=1 and d=2
    if (verbose==0) return(!any(duplicated(unlist(hilf, recursive = FALSE)))) else{
       aus <- !any(duplicated(unlist(hilf, recursive = FALSE)))
       attr(aus, "liste") <- hilf
       return(aus)
      }
  }
  ##################

  ##################
  if (d==2 && t==2){
    ## have to cover pairs of 2-way interactions
    sammel <- array(data=vector(mode="list"),
                    dim=c(choose(m,t), choose(m,t)))
    cases <- nchoosek(m, t)  ## column numbers refer to D
    if (verbose>0)
      rownames(sammel) <- colnames(sammel) <-
      apply(cases, 2, function(obj) paste(obj, collapse="-"))
    anzahlen <- lengths(allrhos)
    ## numbers of instances of the t-way interactions
    ## length of the vector is choose(m,t)
    nintclasses <- length(anzahlen)

    ## lengths of single vectors (for upto=TRUE)
    #vec0s <- mapply(function(obj) 1:obj - start0, levels.no(D))
    ## index number of interaction
    vecs <- mapply(function(obj) 1:obj, anzahlen)

    dset <- nchoosek(nintclasses, d)
    ## add cases with same interaction columns (but different instances)
    ## perhaps do with gtools::combinations and repeats.allowed=TRUE later
    dset <- cbind(dset, rbind(1:nintclasses, 1:nintclasses))
    ## nintclasses is choose(m,t) (corresponding to pairs in cases)
    ##     each column of dset holds two elements of 1:nintclasses
    ##     they can also be the same
    ## nintcombis for d (i.e., if not dupto)
    nintcombis <- ncol(dset)
    for (i in 1:nintcombis){
      cols <- dset[,i]  ## a particular interaction combination (e.g. colums (3,4) with columns (2,5))
      {
        x <- do.call(expand.grid, vecs[cols])  ## all actual variants combined,
                                               ## e.g. {(3,0), (4,1)} with {(2,1),(5,2)}
               ## the selected combinations, elements from interaction indices
               ## has only d=2 columns
        if (cols[1]==cols[2]) {
          ## keep each interaction only once
          ## remove duplicated same interaction, if dupto is FALSE
          if (dupto)
            x <- x[x[,1] <= x[,2],]  ## same covers d=1
          else
            x <- x[x[,1] < x[,2],]
        }
        sammel[min(cols), max(cols)][[1]] <-
          lapply(1:nrow(x), ## every combination
                function(obj){
                  ## obj is the row no of x
                  ## obj2 is the list of rho sets for the tfis of dset[,i],
                  ## j is the list of selections from each list element, as indicated in x
                  ## create list of row number lists
                  hilf <- mapply(function(obj2, j) obj2[[j]],
                                 obj2=allrhos[cols],
                                 j=x[obj,],
                                 SIMPLIFY = FALSE)
                  sort(Reduce(union, hilf))
                  })
        if (verbose>0) names(sammel[min(cols), max(cols)][[1]]) <-
          apply(as.matrix(x), 1, function(obj) paste(obj, collapse=":"))
        ## these names refer to the interaction number
      }
    }
    hilf <- c(sammel)  ## array structure removed
    if (verbose>0) {
      hilf2 <- expand.grid(rownames(sammel), colnames(sammel))
      names(hilf) <- paste(hilf2[,1], hilf2[,2], sep="-")
    }
    hilf <- hilf[sapply(hilf, function(obj) !is.null(obj))]  ## only diagonal and above populated
    ## dupto was already covered by allowing the same interaction in the loop over dset columns
    # if(dupto){
    #   allrhos <- rhos_for_ints(D, t, upto=upto, start0=start0)
    #   ## a list-valued entry for each t-factor combination
    #   ##     length of the list element is the number of possible instances
    #   ##     of interactions for the specific set of t factors
    #   hilf <- c(unlist(allrhos, recursive=FALSE), hilf)
    # }
    if (verbose==0)
    return(!any(duplicated(unlist(hilf, recursive=FALSE))))   ## t=2 and d=2
    else{
      aus <- !any(duplicated(unlist(hilf, recursive=FALSE)))
      attr(aus, "liste") <- hilf
      return(aus)
    }
  }else{"This combination of d and t has not yet been implemented."}
}

is.DA <- function(D, t, d, upto=FALSE, start0=TRUE, dupto=FALSE, verbose=0){
  hilf <- is.LA(D, t, d, upto=upto, start0=start0, dupto=dupto, verbose=1)
  if (!hilf){
    if (verbose>0) return(hilf) else return(FALSE)
  }
  hilf <- unlist(attr(hilf, "liste"), recursive = FALSE)
  ## perhaps have unlisted already as input, or do not unlist here either
  ordnen <- sort(lengths(hilf), index.return=TRUE)
  laengen <- ordnen$x
  hilford <- hilf[ordnen$ix]
  hilford <- lapply(hilford, as.set)
  ## it has to be checked whether any of the sets is a proper subset of any other
  ## there cannot be equal sets because of the LA check
  isDA <- TRUE
  ll <- unique(laengen)
  lluse <- setdiff(ll, max(ll))
  for (lnow in lluse){
    hilfl <- hilford[which(laengen==lnow)]
    hilflonger <- hilford[which(laengen>lnow)]
    for (ik in 1:length(hilfl))
      for (il in 1:length(hilflonger))
        if (set_is_proper_subset(hilfl[[ik]], hilflonger[[il]])){
          isDA <- FALSE
          break
        }
    if (!isDA) break
  }
  if (verbose==0) return(isDA) else{
    attr(isDA, "liste") <- hilf
    return(isDA)
  }
}

metrics <- function(D, t, constr=NULL, start0=TRUE){
  if (!is.null(constr)) stop("constr has not yet been implemented")
  auswert <- is.CA(D, t, verbose=1, start0 = start0)
  lls <- levels.no(D)
  m <- ncol(D)
  hilf <- nchoosek(m,t)
  denoms_diversity <-
    pmin(nrow(D), sapply(1:ncol(hilf), function(obj)
    prod(lls[hilf[,obj]])))
  no_of_tFIs <- 56*choose(m,t)
    ##    sum(denoms_diversity)  ## maximum possible number of different 2fis
    ## given the number of runs and level pattern
  hilf <- attr(auswert, "frequDistr")
  hilf2 <- attr(auswert, "proportions")
  if (names(hilf)[1]=="0"){
    ave.cover <- mean(hilf2) ## average coverage proportion for t-projections
    min.cover <- min(hilf2)
    prop.tot.cover <- mean(hilf2==1)  ## percentage fully-covered t-projections
    cover <- sum(hilf[-1])/sum(hilf)  ## percentage of covered tFIs
    divers <- sum(hilf[-1])/no_of_tFIs #sum(denoms_diversity)
    divers.real <- sum(hilf[-1])/sum(denoms_diversity)
  } else{
    ave.cover <- min.cover <- prop.tot.cover <- cover <- 1
    divers <- sum(hilf)/no_of_tFIs #sum(denoms_diversity)
    divers.real <- sum(hilf)/sum(denoms_diversity)
  }
  to_be_covered <- sum(hilf)
  return(c(cover=cover, to_be_covered=to_be_covered,
           prop.tot.cover=prop.tot.cover, ave.cover=ave.cover,
           min.cover=min.cover,
           projections=choose(m,t),
           divers=divers, no_of_tFIs=no_of_tFIs,
           divers.real=divers.real, no.real=sum(denoms_diversity)))
}

proportions_slower <- function(D, t, verbose=0){
  ## slow because of the overhead in function table
  lls <- levels.no.NA(D)
  m <- ncol(D)
  projs <- nchoosek(m,t)
  nproj <- ncol(projs)
  tots <- sapply(1:nproj, function(obj) prod(lls[projs[,obj]]))
  tabs <- lapply(1:nproj, function(obj) do.call(table,
                        lapply(projs[,obj], function(obj) D[,obj])))
  ncovereds <- sapply(tabs, function(obj) sum(obj>0))
  whichnotcovereds <- lapply(tabs, function(obj) which(obj==0, arr.ind=TRUE))
  total <- sum(ncovereds)/sum(tots)
  proportions <- ncovereds/tots
  ave <- mean(proportions)
  min <- min(proportions)
  simple <- mean(ncovereds==tots)
  if (verbose==1)
    return(list(total=total, ave=ave, min=min, simple=simple,
              tots=tots, ncovereds=ncovereds, proportions=proportions))
  if (verbose==2)
    return(list(total=total, ave=ave, min=min, simple=simple,
                tots=tots, ncovereds=ncovereds,
                proportions=proportions,
                projections=projs,
                tabs=tabs,
                notcovereds=whichnotcovereds))
  ## remaining cases (anything but 1 and 2)
  return(list(total=total, ave=ave, min=min, simple=simple))
}

fasttab <- function(mat,
                    dnn=colnames(mat)){
  ## mat is a matrix whose columns hold integer values
  ## starting with 1
  ## that are to be tabulated
  ## missing values are ignored
  stopifnot(is.matrix(mat))
  stopifnot(all(mat%%1==0, na.rm=TRUE))
  stopifnot(all(mat>0, na.rm=TRUE))
  stopifnot(min(mat, na.rm=TRUE)==1)
  if (ncol(mat)==1) aus <- tabulate(mat)
  else{
    bin <- 0L
    k <- ncol(mat)
    dims <- levels.no.NA(mat)
    pd <- cumprod(c(1,dims))
    dn <- mapply(":", rep(1,k),dims, SIMPLIFY = FALSE)
    for (i in 1:k){
      ## ith column
      a <- mat[,i]
      bin <- bin + pd[i]*(a - 1L)
    }
    if (length(bin)) bin <- bin + 1L
    names(dn) <- dnn
    aus <- array(tabulate(bin, pd[k+1]), dims, dimnames = dn)
    class(aus) <- "table"
  }
  aus
}

LS_to_D3cols <- function(LS, start0=TRUE){
  ## LS a square v x v matrix with entries from 0 to (v-1) or from 1 to v
  ## start0 governs the output: which coding should D3cols have?
  stopifnot(is.matrix(LS))
  v <- nrow(LS)
  stopifnot(v==ncol(LS) & v==length(table(LS)))
  stopifnot(length(unique(table(LS)))==1)
  stopifnot(all(LS %% 1 == 0))  ## integer entries
  stopifnot(min(LS)>0 && max(LS)<=v)
  if (min(LS)==1 && start0) LS <- LS - 1
  if (min(LS)==0 && !start0) LS <- LS + 1
  if (start0) levs <- 0:(v-1) else levs <- 1:v
  cbind(A=rep(levs, each=v), B=rep(levs, v), C=c(LS))
}

D3cols_to_LS <- function(D, start0=NULL){
  ## D a latin square in three column notation (matrix)
  ##   first column taken as column number
  ##   second column taken as row number
  ##   third column taken as entry
  ## start0 logical that governs the output,
  ##   input is coded as 0 to v-1 or 1 to v
  ##   for NULL, same as input
  stopifnot(is.matrix(D))
  stopifnot(is.numeric(D))
  stopifnot(all(D %% 1 == 0))
  levs <- levels.no(D)
  stopifnot(length(unique(levs))==1)
  v <- levs[1]
  stopifnot(min(D)>=0)
  stopifnot(max(D) <= v - 1 + min(D))
  if (is.null(start0)) start0 <- min(D)==0
  if (min(D)==0) D <- D + 1
  stopifnot(all(DoE.base::GWLP(D)==c(1,0,0,v-1)))
  LS <- matrix(NA,v,v)
  for (j in 1:v)
    for (i in 1:v)
      LS[i,j] <- D[D[,1]==j & D[,2]==i ,3]
  if (start0) LS <- LS - 1
  dimnames(LS) <- list(rows=1:v-start0, columns=1:v-start0)
  LS
}

