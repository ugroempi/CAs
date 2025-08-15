#' auxiliary functions for mixed level CAs
#'
#' non-visible functions to be used in functions for MCA creation,
#' and currently also previous drafts that may be deleted
#' key resource: Sherwood (2008) and Groemping (2025)
#'
#' @rdname auxiliary_mixed
#'
#' @usage NULL
#'
#' @details
#' OD: ingredient creation for array expansion
#' one_to_howmany: determines the number of factors a column can be expanded to
#' colexpandOD: expands a column
#' makeNA: makes appropriate levels into flexible (NA) values
#' colexpandODr: repeat colexpandOD
#' find_mix_construction_separate: finds construction details for
#'    a requested setting of levels, for use in MCA2
#' find_mix_construction: unfinished, attempts to be more efficient
#'    but is way more complicated
#' outerseparate: workhorse function of function MCA2

one_to_howmany <- function(N, vfrom, vto, ...){
  ## if the ingoing design is not balanced, use most frequent levels
  ## for making them flexible
  ## assuming near balance, i.e., at most one difference
  perlevel <- ceiling(N/vfrom)
  ## nperlevel says how many levels there are with perlevel entries
  if (!perlevel==N/vfrom) nperlevel <- N%%vfrom else nperlevel <- vfrom
  if (vfrom - vto <= nperlevel)
    Nflex <- perlevel*(vfrom-vto)
  else
    Nflex <- nperlevel*perlevel + (vfrom-vto-nperlevel)*floor(N/vfrom)
  r <- floor(Nflex/(vto*(vto-1)))
  ncol(OD2(vto))^r
}
one_to_howmany_vec <- Vectorize(one_to_howmany, vectorize.args = c("vfrom", "vto"))
# one_to_howmany(25,5,3)

colexpandOD <- function(D, c, lev=2){
  ## D a design whose levels were already reduced to the desired target lev
  ## c a column of that design
  ## internal, for use within a construction
  levs <- levels.no.NA(D)
  stopifnot(levs[c]==lev)
  nas <- which(is.na(D[,c]))
  stopifnot(length(nas)>=lev*(lev-1))
  OD <- OD2(lev)
  new <- do.call(cbind, rep(list(D[,c]), ncol(OD)))
  new[nas[1:(lev*(lev-1))],] <- OD
  if (c==1) aus <- cbind(D[,-1], new)
  if (c==ncol(D)) aus <- cbind(D[,-ncol(D)], new)
  if (1 < c && ncol(D)>c) aus <- cbind(D[,c(1:(c-1), (c+1):ncol(D))], new)
  aus
}

makeNA <- function(D, levs){
  # D a matrix
  # levs a vector of levels
  D[which(D %in% levs),] <- NA
  D
}

colexpandODr <- function(D, c, lev=2, ncolmax=Inf){
  ## repeat the process as long as possible
  mat <- D[,c, drop=FALSE]
  r <- sum(is.na(mat)) %/% (lev*(lev-1))
  for (i in 1:r){
    for (cc in ncol(mat):1)
      mat <- colexpandOD(mat, cc, lev)
    if (ncol(mat) >= ncolmax) break
  }
  ncolmax <- min(ncol(mat), ncolmax)
  mat[,1:ncolmax, drop=FALSE]
}

find_mix_construction_separate <- function(nlevels, ...){
  ## nlevels is a data.frame with a column levels and a column frequency
  ## strength 2 only
  if (!all(sort(nlevels$levels, decreasing = TRUE)==nlevels$levels))
    stop("nlevels must be sorted in decreasing order")
  vmax <- nlevels$levels[1]
  ## prepare matrix for output object
  ## which is the same for all cases in this crude separate implementation
  assigncols <- diag(nlevels$frequency)
  dimnames(assigncols) <- list(donor=nlevels$levels, receiver=nlevels$levels)

  ## estimate N (might be improvable by allowing overlap between levels)
  ## separate column for each number of levels
  N <- bestN(2, nlevels$frequency[1] + length(nlevels$levels)-1, vmax)
  ## for this N, what's needed?
  percolumn <- sapply(nlevels$levels, function(obj) one_to_howmany(N, vmax, obj))
  columnsneeded_separate <- ceiling(nlevels$frequency/percolumn)
  names(columnsneeded_separate) <- nlevels$levels
  ## naively requested columns insufficient ?
  if (sum(columnsneeded_separate) > nlevels$frequency[1] + length(nlevels$levels)-1){
    ## case where the run size is sufficient also for the number of columns that is needed
    if (bestN(2, sum(columnsneeded_separate), vmax)==N)
      return(list(setting=nlevels, k=sum(columnsneeded_separate), N=N, v=vmax,
                  expand=list(ingoing_needed=columnsneeded_separate, assign=assigncols)))
    ## N is not enough, increase
    Nold <- N
    N <- bestN(2, sum(columnsneeded_separate), vmax)
    # print("eins")
    # print(N)
    percolumn_new <- sapply(nlevels$levels, function(obj) one_to_howmany(N, vmax, obj))
    columnsneeded_separate_new <- ceiling(nlevels$frequency/percolumn_new)
    names(columnsneeded_separate_new) <- nlevels$levels
    ## decreased number of needed columns versus prior N
    ## more columns available than needed when reducing the requested k in bestN
    if (sum(columnsneeded_separate_new) < sum(columnsneeded_separate) &&
        sum(columnsneeded_separate_new) > nlevels$frequency[1] + length(nlevels$levels)-1){
      ## maybe columnsneeded_separate_new is too large, implying a too large N
      ## therefore another adjustment attempt
      Noldlarge <- N
      N <- bestN(2, sum(columnsneeded_separate_new), vmax)
      # print("zwei")
      # print(N)
      percolumn_new2 <- sapply(nlevels$levels, function(obj) one_to_howmany(N, vmax, obj))
      columnsneeded_separate_new2 <- ceiling(nlevels$frequency/percolumn_new2)
      names(columnsneeded_separate_new2) <- nlevels$levels
      if (sum(columnsneeded_separate_new2) >= nlevels$frequency[1] + length(nlevels$levels)-1)
        return(list(setting=nlevels, k=sum(columnsneeded_separate_new2), N=N, v=vmax,
                    expand=list(ingoing_needed=columnsneeded_separate_new2, assign=assigncols)))
      else
        return(list(k=sum(columnsneeded_separate_new), N=Noldlarge, v=vmax,
                    expand=list(ingoing_needed=columnsneeded_separate_new, assign=assigncols)))
    }else
      return(list(setting=nlevels, k=sum(columnsneeded_separate), N=N, v=vmax,
                  expand=list(ingoing_needed=columnsneeded_separate,
                              assign=assigncols)))
  }else
    return(list(setting=nlevels, k=sum(columnsneeded_separate), N=N, v=vmax,
                expand=list(ingoing_needed=columnsneeded_separate,
                            assign=assigncols)))
}



## find_mix_construction_separate(data.frame(levels=c(8,4,3,2), frequency=c(1,5,7,12)))
## find_mix_construction_separate(data.frame(levels=c(6,4,3,2), frequency=c(2,8,30,70)))

find_mix_construction <- function(nlevels, ...){
  ## nlevels is a data.frame with a column levels and a column frequency
  ## strength 2 only

  ## the expand element of the output object is a list whose elements are
  ##   - a vector with numbers of columns, like in find_mix_construction_separate
  ##   - and a matrix with numbers of resulting columns attributed to ingoing columns
  ##     from different numbers of levels: row = donor levels, col = receiver levels
  if (!all(sort(nlevels$levels, decreasing = TRUE)==nlevels$levels))
    stop("nlevels must be sorted in decreasing order")
  vmax <- nlevels$levels[1]
  names(nlevels$frequency) <- nlevels$levels

  hilf <- find_mix_construction_separate(nlevels, ...)

  ## goal: reduce k
  for (j in 1:(hilf$k-2)){
    N <- bestN(2, hilf$k - j, vmax)
    if (N < hilf$N) break
  }
  know <- hilf$k - j ## new k
print(know) ## 4, NIST 43 runs

  ## for this N, what's needed?

  # percolumn_mat <- outer(nlevels$levels, nlevels$levels, function(X, Y) one_to_howmany_vec(N, X, Y))
  # rownames(percolumn_mat) <- colnames(percolumn_mat) <- nlevels$levels
  # percolumn <- sapply(nlevels$levels, function(obj) one_to_howmany(N, vmax, obj))

  remainingNAs <- rep(list(0),length(nlevels$levels))
  names(remainingNAs) <- nlevels$levels
  columns_needed <- nlevels$frequency[1]
  remainingNAs[[1]] <- rep(0, columns_needed)
  assigncols <- matrix(0, nrow=length(nlevels$levels), ncol=length(nlevels$levels),
                       dimnames=list(donor=nlevels$levels, receiver=nlevels$levels))
  assigncols[1,1] <- nlevels$frequency[1] ## the largest number of levels is actually needed
  ndonor_cols <- rep(0, length(nlevels$levels))
  ndonor_cols[1] <- columns_needed  ## for vmax, we need an entire column in each case

  for (j in 2:length(nlevels$levels)){
    ## use prior columns as far as possible
    ## update remainingNAs for prior levels entries
    ## create remainingNAs for current levels entry
    ##    (unchanged at 0, if it is created exclusively from prior columns)
    levnow <- nlevels$levels[j]
    ncolODnow <- ncol(OD2(levnow))
    needednow <- nlevels$frequency[j]
    ## how much from leftover NA from previous (=larger) levels?
    remainingNAsprior <- remainingNAs[1:(j-1)]
    ndonor_prior <- ndonor_cols[1:(j-1)]
    ## this is for one donor column of the respective sort
    ## must be a list of numeric vectors, because remainingNAsprior is a list of numeric vector
    potentialColumnsPrior <- lapply(remainingNAsprior,
                                    function(obj) ifelse(levnow*(levnow - 1) > obj,
                                    0, ## no column possible
                                    ncolODnow^floor(obj/(levnow*(levnow - 1)))))
# print(potentialColumnsPrior)
# print(ndonor_prior)
    #stopifnot(all(lengths(potentialColumnsPrior)==ndonor_prior)) ## sanity check

    ## how many NA from a new levnow column before expanding
    initialNAnow <- ceiling(N*(vmax-levnow)/vmax) ## cautious, works with wise choice of levels
                       ## *very* wise choice may occasionally have more NA
    ## how many levnow-level columns can one such column be expanded to?
    potentialColumnsNow <- ncolODnow^floor(initialNAnow/(levnow*(levnow - 1)))
    check <- (needednow - sum(unlist(potentialColumnsPrior))) %/% potentialColumnsNow    ## minimum number of required donor potential columns
    addcheck <- (needednow - sum(unlist(potentialColumnsPrior))) %% potentialColumnsNow
    if (addcheck > 0) check <- check + 1
    ## can all the check required vmax-level columns create enough columns ?
    if (check * potentialColumnsNow >= needednow) useprior <- FALSE else useprior <- TRUE

    ## fill current entries
    ndonor_cols[j] <- check ## may even be zero
    if (sum(ndonor_cols) > know){
      message("use of prior columns does not improve the outcome")
      return(hilf)
    }
    if (useprior){
          assigncols[j,j] <- check * potentialColumnsNow
          needednow <- needednow - check * potentialColumnsNow
          if (check > 0)
            remainingNAs[[j]] <- rep(initialNAnow %% (levnow*(levnow-1)), check * potentialColumnsNow)
          else
            remainingNAs[[j]] <- numeric(0)  ## maybe different coding of this case better?
    }else{## no prior columns used
          ## last column may keep more NAs than first columns
          assigncols[j,j] <- needednow
          remainingNAs[[j]] <- c(rep(initialNAnow %% (levnow*(levnow-1)),
                                         (check-1) * potentialColumnsNow),
                                rep(initialNAnow - levnow*(levnow-1) *
                                  ceiling((needednow %% potentialColumnsNow)/ncolODnow),
                                  needednow %% potentialColumnsNow))
        next  ## move to next j, as no useprior is needed
    }
    ## now we treat prior use
    ## regardless whether there are some new columns or not

        ## allocate needednow columns to donor rows, assigncols[j,j] remains 0 or was already treated
        for (i in (j-1):1){
          ## proceed from smaller to larger levels
          ncolnow <- ndonor_prior[i]              ## number of such columns
          if (ncolnow==0) next
          potcolnow <- potentialColumnsPrior[[i]]   ## vector of entries
          if (sum(potcolnow) <= needednow){
            ## use them all (and likely some more)
            assigncols[i,j] <- assigncols[i,j] + potcolnow
            remainingNA[[i]] <- remainingNA[[i]] - potcolnow/ncolODnow*levnow*(levnow-1)
            if (sum(potcolnow)==needednow) break ## no further need
            needednow <- needednow - sum(potcolnow)
          }else{
            ## last column keeps more NA values
            nneeded <- min(which(cumsum(potcolnow)>=needednow))
            if (nneeded == 1) neededlast <- needednow else
            neededlast <- needednow - sum(potcolnow[-(nneeded:ncol)])
            stopifnot(neededlast < potcolnow[nneeded])
            assigncols[i,j] <- assigncols[i,j] + needednow
            if (nneeded > 1)
            remainingNA[[i]][1:(nneeded-1)] <- remainingNA[[i]][1:(nneeded-1)] -
              potcolnow[1:(nneeded-1)]/ncolODnow*levnow*(levnow-1)
            ## last (or only) needed column
            remainingNA[[i]][nneeded] <- remainingNA[[i]][nneeded] -
              ceiling(neededlast/ncolODnow)*levnow*(levnow-1)
            ## remaining values do not need to be changed
          }
        }
      }  ## end of j loop

  list(setting=nlevels,
       k = sum(ndonor_cols), N=N, v=vmax,
       expand=list(ingoing_needed=ndonor_cols, assign=assigncols))
}

outerseparate <- function(D, nlevels, ...){
  ## nlevels a data.frame with columns levels and frequency
  ## does this work with a CA D, or does it require an OA or
  ##       at least a balanced CA ?
  ## ... not used
  #hilf <- find_mix_construction_separate(nlevels, ...)

  if (any(is.na(D))) stop("D must not have NA values")
  N <- nrow(D); k <- ncol(D); v <- levels.no(D)
  # levels are ordered such that latest levels have largest frequencies
  # which leads to largest possible numbers of NA values in each step
  for (i in 1:ncol(D)){
    tab <- table(D[,i])
    vi <- length(tab)
    if (length(unique(tab)) > 1) D <- permvals(D, i, 0:(vi-1), order((0:(vi-1))[order(tab)])-1)
  }
    #D <- L72.2.8.3.8.4.1.6.5[,18:22] ## eval(parse(text=labelToCode(names(hilf$N), 2, k, v)))
  aus <- matrix(NA, N, sum(nlevels$frequency))
  ## first columns have no reductions
#  aus[,1:nlevels$frequency[1]] <- D[, 1:nlevels$frequency[1]]
  donor_done <- 0         # nlevels$frequency[1]
  targetcols_done <- 0    # donor_done ## for this v only
  nsetting <- nrow(nlevels)

  for (j in 1:nsetting){
    ## columns with v levels skipped
    #ndonornow <-  # hilf$expand$ingoing_needed[j]
    levnow <- nlevels$levels[j]
    ## simple for separate: all from these donor columns
    tocreatenow <- nlevels$frequency[j]
    creatednow <- 0
    while (creatednow < tocreatenow){
      matnow <- matrix(D[,donor_done+1], N, 1)
      ## sort levels so that the most frequent are last
      # levord <- as.numeric(names(sort(table(matnow))))
      # matnow <- permvals(matnow, 1, 0:(v-1), levord)
      matnow[matnow > levnow - 1] <- NA
      if (sum(is.na(matnow)) >= levnow*(levnow-1))
          hilf_expanded <- colexpandODr(matnow, 1, lev=levnow,
                                        ncolmax = tocreatenow)
      else hilf_expanded <- matnow
      donor_done <- donor_done + 1
      if (tocreatenow >= ncol(hilf_expanded)){
        aus[,(targetcols_done+1):(targetcols_done+ncol(hilf_expanded))] <- hilf_expanded
        targetcols_done <- targetcols_done + ncol(hilf_expanded)
      }
      else{
        aus[,(targetcols_done+1):(targetcols_done+tocreatenow)] <- hilf_expanded[,1:tocreatenow]
        targetcols_done <- targetcols_done + tocreatenow
      }
      tocreatenow <- max(0, tocreatenow - ncol(hilf_expanded))
    }## end of while loop
  }## end of j loop
  aus
}
