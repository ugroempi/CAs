#' Functions for checking the Cyclotomy conditions
#'
#' Conditions for constructions 1 to 4 and their modifications
#' of Colbourn (2010) can be checked.
#' The checks take a long time for large cases.
#'
#' @rdname checkCyclotomyConditions
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @aliases checkcond1
#' @aliases checkcond2
#' @aliases checkcond3
#' @aliases checkcond4
#' @aliases checkcond3a
#' @aliases checkcond3b
#' @aliases checkcond4a
#' @aliases checkcond4b
#'
#' @usage checkcond1(t, v, q, iter=FALSE, progress=FALSE)
#' @usage checkcond2(t, v, q, iter=FALSE, progress=FALSE)
#' @usage checkcond3(t, v, q, iter=FALSE, progress=FALSE)
#' @usage checkcond4(t, v, q, iter=FALSE, progress=FALSE)
#' @usage checkcond3a(t, v, q, iter=FALSE, progress=FALSE, type3=TRUE)
#' @usage checkcond3b(t, v, q, iter=FALSE, progress=FALSE, type3=TRUE)
#' @usage checkcond4a(t, v, q, iter=FALSE, progress=FALSE, type4=TRUE)
#' @usage checkcond4b(t, v, q, iter=FALSE, progress=FALSE, type4=TRUE)
#'
#' @param t CA strength to be checked
#' @param v number of levels
#' @param q prime or prime power
#' @param iter logical (default FALSE): if FALSE, all combinations are calculated at once (slightly faster); otherwise, combinations are iterated through one at a time (works for large cases where memory is not sufficient for the other strategy).
#' @param progress logical (default FALSE): if TRUE, show a progress bar.
#' For the a and b types, two progress bars will be shown (second one progresses much faster than first one).
#' @param type3 Should the type 3 check be run? Set to FALSE, if it is already known that the type 3 condition is fulfilled.
#' @param type4 Should the type 4 check be run? Set to FALSE, if it is already known that the type 4 condition is fulfilled.
#'
#' @returns a logical indicating whether the condition is fulfilled.
#'
#' @section Details: \code{checkcond1} to \code{checkcond4} do
#' the basic checks. Checks for constructions 3a and 4a consist of
#' checks for constructions 3 and 4, combined with checks for construction 1
#' for strength t-1,\cr
#' and checks for constructions 3b and 4b consist of checks for constructions 3 and 4,
#' combined with checks for condition 2 for strength t-1.
#'
#' @examples
#' checkcond1(3,2,11)  ## FALSE
#' checkcond2(3,2,11, progress=TRUE)  ## TRUE, example for progress bar
#' ## if construction 4 with strength 4 works,
#' ## then this setting will also fulfill construction 4b
#' checkcond4b(4,2,11, progress=TRUE)
#' N_cyc(11,2,"4b")
#' k_cyc(11,2,"4b")
#' ## this array is the currently best-known
#' eCAN(4,12,2)
#'

#' @export
checkcond1 <- function(t,v,q,iter=FALSE, progress=FALSE){
  ## the implementation with nchoosek does not work
  ## for large values of q and t (e.g., 10007 and 5)
  gf <- lhs::create_galois_field(q)
xstart <- cycvec(v,q,gf=gf)

## every t-subset of GF(q) (sets in columns)
if (!iter)
  tsets <- nchoosek(q, t) ## these start with 1
## every t-combination with replacement in rows
vsets <- as.matrix(do.call(expand.grid,
                           lapply(1:t, function(obj) 0:(v-1))))

if (progress) pb <- txtProgressBar(0, choose(q,t), char=":")

i <- 0
nmax <- choose(q,t)
while(i < nmax){
    i <- i + 1
    if (progress) setTxtProgressBar(pb, i)
    if (iter){
      if (i==1) pick <- 1:t else pick <- gen.next.cbn(pick,q)
    } else pick <- tsets[,i]
    xnow <- xstart[pick]
    for (j in 1:nrow(vsets)){
      fullfillnow <- FALSE ## start for this ij
      vnow <- vsets[j,]
      ## fulfilled with a==0
      if (all(vnow==xnow)) fullfillnow <- TRUE
      else{
        ## check further values of a
        for (a in 1:(q-1)){
          xnowa <- xstart[SOAs:::gf_sum(pick-1,rep(a,t),gf)+1]
          if (all(vnow==xnowa)) fullfillnow <- TRUE
          ## leave a loop after success
          if (fullfillnow) break
        }
      }
      ## go to next j for same i
      if (fullfillnow) next else return(FALSE)
    }
  }

if (progress) close(pb)

## if survived all loop steps, return TRUE
return(TRUE)
}

## needs to become iterative by argument, like checkcond1
## not yet implemented in spite of argument!

#' @export
checkcond2 <- function(t,v,q, iter=FALSE, progress=FALSE){
  ## same as checkcond1, only without the constant non-zero rows of vsets
  gf <- lhs::create_galois_field(q)
  xstart <- cycvec(v,q,gf=gf)
  ## every t-subset of GF(q) (sets in columns)
  if (!iter)
    tsets <- nchoosek(q, t) ## these start with 1
  ## every t-combination with replacement in rows
  vsets <- as.matrix(do.call(expand.grid,
                    lapply(1:t, function(obj) 0:(v-1))))
  ## constant rows with all elements positive
  rem <- which(apply(vsets,1, function(obj) all(obj>0) &&
                       length(table(obj))==1))
  vsets <- vsets[-rem,]

  if (progress) pb <- txtProgressBar(0, choose(q,t), char=":")

  i <- 0
  nmax <- choose(q,t)
  while(i < nmax){
      i <- i + 1
      if (progress) setTxtProgressBar(pb, i)
      if (iter){
        if (i==1) pick <- 1:t else pick <- gen.next.cbn(pick,q)
      } else pick <- tsets[,i]
    xnow <- xstart[pick]
    for (j in 1:nrow(vsets)){
      fullfillnow <- FALSE ## start for this ij
      vnow <- vsets[j,]
      ## fulfilled with a==0
      if (all(vnow==xnow)) fullfillnow <- TRUE
      else{
        ## check further values of a
        for (a in 1:(q-1)){
          xnowa <- xstart[SOAs:::gf_sum(pick-1,rep(a,t),gf)+1]
          if (all(vnow==xnowa)) fullfillnow <- TRUE
          ## leave a loop after success
          if (fullfillnow) break
        }
      }
      ## go to next j for same i
      if (fullfillnow) next else return(FALSE)
    }
  }
    ## if survived all loop steps, return TRUE
    if (progress) close(pb)
    return(TRUE)
}

## needs to become iterative by argument, like checkcond1
## not yet implemented in spite of argument!

#' @export
checkcond3 <- function(t,v,q,iter=FALSE, progress=FALSE){
  ## 3a cases must pass this check
  ## i.e., checkcond3(4,2,23) must yield TRUE
  gf <- lhs::create_galois_field(q)
  xstart <- cycvec(v,q,gf=gf)
  ## every t-subset of GF(q) (sets in columns)
  if (!iter)
    tsets <- nchoosek(q, t) ## these start with 1
  ## every t-combination with replacement in rows
  vsets <- as.matrix(do.call(expand.grid, lapply(1:t, function(obj) 0:(v-1))))

  if (progress) pb <- txtProgressBar(0, choose(q,t), char=":")

  i <- 0
  nmax <- choose(q,t)
  while(i < nmax){
    i <- i + 1
    if (progress) setTxtProgressBar(pb, i)
    if (iter){
      if (i==1) pick <- 1:t else pick <- gen.next.cbn(pick,q)
    } else pick <- tsets[,i]
    xnow <- xstart[pick]
    for (j in 1:nrow(vsets)){
      fullfillnow <- FALSE ## start for this ij
      vnow <- vsets[j,]
      ## fulfilled with a==0 and r==0
      if (all(vnow==xnow)) fullfillnow <- TRUE
      else{
        ## check further values of r for a==0
        for (r in 1:(v-1)){
          vnowr <- (vnow+r)%%v
          if (all(vnowr==xnow)) fullfillnow <- TRUE
          ## leave r loop after success
          if (fullfillnow) break
        }
        ## check further values of a
        if (!fullfillnow){
        for (a in 1:(q-1)){
          xnowa <- xstart[SOAs:::gf_sum(pick-1,rep(a,t),gf)+1]
          ## check for r==0
          if (all(vnow==xnowa)) fullfillnow <- TRUE
          else{
            ## check further values of r
            for (r in 1:(v-1)){
              vnowr <- (vnow+r)%%v
              if (all(vnowr==xnowa)) fullfillnow <- TRUE
              ## leave r loop after success
              if (fullfillnow) break
              }
          }
          ## leave a loop after success
          if (fullfillnow) break
        }
      }}
      ## go to next j for same i
      if (fullfillnow) next else return(FALSE)
    }
  }
  if (progress) close(pb)
  ## if survived all loop steps, return TRUE
  return(TRUE)
}
## checkcond3(4,2,23) ## TRUE expected

#' @export
checkcond3a <- function(t,v,q,iter=FALSE, progress=FALSE, type3=TRUE){
  if (type3){
    if (!checkcond3(t,v,q,iter=iter, progress=progress)) return(FALSE)
    else(print("type 3 successful, continue with check of type 1 for t-1"))
  }
  checkcond1(t-1,v,q,progress=progress)
}
## checkcond3a(4,2,23) ## TRUE expected

#' @export
checkcond3b <- function(t,v,q,iter=FALSE, progress=FALSE, type3=TRUE){
  if (type3){
  if (!checkcond3(t,v,q,iter=iter, progress=progress)) return(FALSE)
    else(print("type 3 successful, continue with check of type 2 for t-1"))
  }
  checkcond2(t-1,v,q,iter=iter, progress=progress)
}
## 1 cases are also TRUE here
## checkcond3b(3,2,19)  ## TRUE expected
## checkcond3b(3,4,41)  ## paper claims TRUE, but found FALSE in practice
##                         hence, FALSE seems correct
## checkcond3b(5,2,103) ## cyc(104,5,2,103,type="3b")
##     did not complete but ran a long time (more than 24h) without returning FALSE


## needs to become iterative by argument, like checkcond1
#' @export
checkcond4 <- function(t,v,q,iter=FALSE, progress=FALSE){
  ## that is 3, but without the constant rows in vsets
  ##              (regardless of zero or non-zero)
  ## 4a cases must pass this check
  ##  checkcond4(3,3,19) should yield TRUE
  ##  checkcond4(4,2,19) should yield TRUE
  ##
  gf <- lhs::create_galois_field(q)
  xstart <- cycvec(v,q,gf=gf)
  ## every t-subset of GF(q) (sets in columns)
  if (!iter)
    tsets <- nchoosek(q, t) ## these start with 1
  ## every t-combination with replacement in rows
  vsets <- as.matrix(do.call(expand.grid, lapply(1:t, function(obj) 0:(v-1))))
  rem <- which(apply(vsets,1, function(obj) length(table(obj))==1))
  vsets <- vsets[-rem,]

  if (progress) pb <- txtProgressBar(0, choose(q,t), char=":")

  i <- 0
  nmax <- choose(q,t)
  while(i < nmax){
    i <- i + 1
    if (progress) setTxtProgressBar(pb, i)

    if (iter){
      if (i==1) pick <- 1:t else pick <- gen.next.cbn(pick,q)
    } else pick <- tsets[,i]
    xnow <- xstart[pick]
    for (j in 1:nrow(vsets)){
      fullfillnow <- FALSE ## start for this ij
      vnow <- vsets[j,]
      ## fulfilled with a==0 and r==0
      if (all(vnow==xnow)) fullfillnow <- TRUE
      else{
        ## check further values of r for a==0
        for (r in 1:(v-1)){
          vnowr <- (vnow+r)%%v
          if (all(vnowr==xnow)) fullfillnow <- TRUE
          ## leave r loop after success
          if (fullfillnow) break
        }
        ## check further values of a
        if (!fullfillnow){
          for (a in 1:(q-1)){
            xnowa <- xstart[SOAs:::gf_sum(pick-1,rep(a,t),gf)+1]
            ## check for r==0
            if (all(vnow==xnowa)) fullfillnow <- TRUE
            else{
              ## check further values of r
              for (r in 1:(v-1)){
                vnowr <- (vnow+r)%%v
                if (all(vnowr==xnowa)) fullfillnow <- TRUE
                ## leave r loop after success
                if (fullfillnow) break
              }
            }
            ## leave a loop after success
            if (fullfillnow) break
          }
        }}
      ## go to next j for same i
      if (fullfillnow) next else return(FALSE)
    }
  }
  if (progress) close(pb)

  ## if survived all loop steps, return TRUE
  return(TRUE)
}

#' @export
checkcond4a <- function(t,v,q,iter=FALSE, progress=FALSE, type4=TRUE){
  ##  checkcond4a(3,3,19) should yield TRUE
  ##  checkcond4a(4,2,19) should yield TRUE
  ##  checkcond4a(3,4,37) should yield TRUE
  if (type4){
    if (!checkcond4(t,v,q, iter=iter, progress=progress)) return(FALSE)
    else(print("type 4 successful, continue with check of type 1 for t-1"))
  }
  checkcond1(t-1,v,q,iter=iter, progress=progress)
}

#' @export
checkcond4b <- function(t,v,q,iter=FALSE, progress=FALSE, type4=TRUE){
  if (type4){
    if (!checkcond4(t,v,q, iter=iter, progress=progress)) return(FALSE)
    else(print("type 4 successful, continue with check of type 2 for t-1"))
  }
  checkcond2(t-1,v,q,iter=iter, progress=progress)
}
