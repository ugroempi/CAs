#' coverage and diagnostic plot functions
#'
#' functions for calculating covering proportions for a given strength,
#' and for visualizing those
#'
#' @rdname coverage
#' @aliases coverage
#' @aliases coverage_iter
#' @aliases prepcoverplot
#' @aliases coverplot
#'
#' @usage coverage(D, t, isInteger=TRUE, verbose=0, start0=TRUE, parallel=1)
#' @usage coverage_iter(D, t, isInteger=TRUE, start0=TRUE,
#'        progress=FALSE, abortnot1=FALSE, start.proj=1)
#' @usage coverplot(D, t, isInteger=TRUE, start0=TRUE,
#'     type="projections",
#'     main=NULL, xlab=NULL, ylab=NULL,
#'     col=rgb(0,0,1,0.5), las=1, plot=TRUE, ...)
#'
#' @param D array
#' @param t small positive integer, interaction strength to be considered
#' @param isInteger logical (default \code{TRUE}): are the elements of D all
#'                  consecutive integers starting with 0 or 1?
#' @param start0 logical (default \code{TRUE}): do the integers start with 0 ? (irrelevant for \code{!isInteger})
#' @param verbose 0 (default) returns only four percentages, 1 and 2 return more detail.
#' Anything apart from 1 and 2 behaves like 0.
#' @param parallel number of threads to be used (default: 1); if increased, package parallel must be available
#' @param progress logical (default \code{FALSE}): should a progress bar be printed ?
#' @param abortnot1 logical (default \code{FALSE}): should the function stop at the first projection for which
#' strength \code{t} coverage is not perfect ?
#' @param start.proj integer, the number of the projection to start with in lexicographical
#' order; ignored for \code{abortnot1=FALSE}; for \code{abortnot1=TRUE}, assumes that all
#' projections up to \code{start.proj - 1} were already successfully checked
#' @param type \code{"projections"} (default) or \code{"tuples"}
#' @param main title of plot
#' @param xlab,ylab axis titles
#' @param las orientation of tickmark labels
#' @param col color
#' @param plot logical (default \code{TRUE}): generate the plot?
#' @param ... further arguments to polygon
#'
#' @return \code{coverage} returns an object of class coverage, which is a
#' list with elements \code{total}, \code{ave}, \code{min} and \code{simple},
#' with increasing detail added for verbose 1 or 2 (function \code{coverage}) or
#' where a run was aborted because of \code{abortnot1=TRUE} (function \code{coverage_iter}).
#'
#' The \code{coverplot} functions return a list of \code{x}, \code{y}, \code{t}
#' and \code{proportions} (from call with \code{verbose=1}).
#'
#' @section Details:
#' The proportions calculated are the proportion of covered
#' \code{t}-tuples (\code{total}, taking the tuples as units), the average proportion of covered
#' \code{t}-tuples taken over the \code{t}-factor projections (\code{ave}),
#' the \code{min}imum proportion of \code{t}-tuples covered taken over the \code{t}-factor
#' projections, and the proportion of \code{t}-factor projections that are
#' completely covered (\code{simple}).
#'
#' For function \code{coverage}, if verbosity is requested, there are additional list elements:
#'
#' for \code{verbose=1} the number of covered tuples for each projection (\code{ncovereds}),
#' the proportion this number represents (\code{proportions}),
#'
#' for \code{verbose=2} the frequency tables of tuples (\code{tabs}),\cr
#' a list of identification indices of uncovered tuples for each projection (\code{uncovereds})\cr
#' and the \code{projections} as a \code{t x choose(k, t)}-matrix that can be used as the legend
#' for the columns of \code{tabs} for the elements of \code{uncovered}.
#'
#' Function \code{coverage_iter} calculates the same quantities as function \code{coverage};
#' instead of creating all projections at once, it creates them iteratively, which is slower
#' but doable for large cases for which function \code{coverage} fails for memory reasons.
#' This iterative version is, at least at present, not parallelized. It allows a progress bar
#' (which slightly deteriorates run time, e.g., on my Windows 10 machine, about 7 milliseconds
#' per 1000 projections). It also allows to request abortion on encountering the first violation
#' of strength \code{t} coverage.\cr
#' In case of aborting, the function returns \code{NA} values for
#' all four metrics, and as a fifth list element, the message at which projection the violation
#' occurred; \code{coverage_iter} does not have a \code{verbose} argument, but using the
#' \code{verbose} argument on printing an object that resulted from an aborted run prints
#' that message.\cr
#' For \code{abortnot1=TRUE}, function runs that did not abort up to a particular
#' iteration can be restarted at the subsequent iteration. The actual restart point will be
#' the first projection starting with the same first element as the requested \code{start.proj}
#' e.g., for \code{t=3} and a matrix \code{D} with 6 columns, the first ten 3-factor projections
#' start with 1, followed by six projections that start with 2, three projections that start with 3
#' and a single projection that starts with 4 (10 + 6 + 3 + 1 = 20). \code{start.proj=12} will start
#' with the projection 2:3:4, i.e., the 11th projection,
#' which is the first that starts with 2. For projections that start with 1, the second factor is
#' also relevant for the start position: the first four projections have 2 in the second position,
#' then there are three with 3 in the second position, then two with 4, and finally one with 5.
#' Thus, for \code{start.proj=7}, \code{start.proj} will be reduced to 5, and the projection
#' is 1:3:4 (the first that starts with 1:3).
#'
#' Per default, the print method for class \code{coverage} prints the first four
#' elements only.
#'
#' Function \code{coverplot}, with \code{type=projections}, plots the
#' percentage of projections on the horizontal axis and
#' the coverage on the vertical axis (visualizing the average taken over projections);
#' with \code{type=tuples}, the percentage of tuples (i.e., projections as weighted by tuples)
#' is shown on the horizontal axis, visualizing the total percentage covered.
#' Sorting of the combinations on the horizontal axis is from best to worst coverage.
#'
#' @section Warning:
#' For function \code{coverage_iter}, in case of starting at projection \code{start.now},
#' the user is responsible to ensure that all prior runs have also been checked; the calculated
#' values assume that the prior runs showed perfect coverage.
#'
#' @examples
#' ## fixed levels
#' D <- KSK(k=10)
#' coverage(D, 2)
#' coverage(D, 3)
#' ## mixed levels, total and ave are different
#' plan <- cbind(rep(0:2, each=3), c(rep(0:1, 4), 1), c(rep(0:1, each=4),0))
#' coverage(plan, 2)
#' ## verbose content not printed per default
#' coverage(plan, 2, verbose=1)
#' print(coverage(plan, 2, verbose=1), digits=4, verbose=TRUE)
#' print(coverage(plan, 2, verbose=2), digits=4, verbose=TRUE)
#' ## two of three projections fully covered,
#' ## one covered 5/6 only
#'
#' coverplot(plan, 2, type="projections")
#' ## the not fully covered projection involves the 3-level
#' ## column and thus has 6/(6+6+4)=37.5 pct of the tuples instead of 1/3
#' coverplot(plan, 2, type="tuples")
#'
#' ## the Hadamard based 16 run design for 14 columns
#' ##    that is available from FrF2::pb;
#' ##    the design has levels -1 and 1 and columns are factors
#' ##    therefore, isInteger=FALSE has to be used
#' if (require(FrF2)){
#'    CA3.2tothe14 <- FrF2::pb(16, 14)
#'    print(coverage(CA3.2tothe14, 3, isInteger=FALSE))
#'    }
#'


#' @export
coverage <- function(D, t, isInteger=TRUE,
                        verbose=0, start0=TRUE, parallel=1){
  ## made faster by custom table function fasttab
  ## isInteger: are the values all from 0 to v-1 or from 1 to v?
  ## start0: 0 to v-1? (relevant for isInteger==TRUE only)
      ## parallel: integer number of threads
      ## for t=4 with pb(48), 4 threads are very slightly slower than 7 threads and
      ##       take about half the time of 1 thread
      ## for t=5, the advantage of 7 over 4 threads is slightly larger
      ##      62.46        3.00      174.88 from 7 threads
      ##      48.84        2.32      190.78 from 4 threads (machine has 4 physical threads)
      ##     308.45        8.58      317.34 from a single thread
  ## loop variant was as fast or slightly faster for 1 thread,
  ##    but much slower with foreach and doParallel
  ## thus, large cases are likely not doable
  stopifnot(is.matrix(D) || is.data.frame(D))
  if (is.data.frame(D) && isInteger) D <- as.matrix(D)
  if (is.matrix(D) && !isInteger) D <- as.data.frame(D)
  if (is.matrix(D)) if (!is.numeric(D))
    stop("isInteger is TRUE, but D has non-numeric content")
  if (is.matrix(D) && start0){
    if (!all(apply(D, 2, min)==0))
      stop("start0 is TRUE, but columns do not start at 0")
    D <- D + 1
  }
  if (is.matrix(D) && !start0){
    if (!all(apply(D, 2, min)==1))
      stop("start0 is FALSE, but columns do not start at 1")
  }
  m <- ncol(D)
  if (is.data.frame(D)){
    ## create integers via factors
    ## (also treats matrices for isInteger=FALSE)
    for (i in 1:m){
      D[[i]] <- as.integer(factor(D[[i]]))
    }
    D <- as.matrix(D)
  }
  ## now D is a matrix with values starting at 1
  ## (provided that the integer versions are adequate)
  lls <- levels.no.NA(D)  ## does not count NA values
  projs <- nchoosek(m,t)
  # nproj <- ncol(projs)

  tots <- apply(projs, 2, function(obj) prod(lls[obj]))
  if (parallel==1)
    tabs <- apply(projs, 2, function(obj) fasttab(D[,obj, drop=FALSE]), simplify = TRUE)
  else{
    stopifnot(requireNamespace("parallel"))
    stopifnot(parallel <= parallel::detectCores())
    mycl <- parallel::makePSOCKcluster(parallel)
    parallel::clusterExport(mycl, c("projs", "D"), envir=environment())
    parallel::clusterExport(mycl, c("fasttab"), envir=environment(CAs:::fasttab))
#    tabs <- parallel::parLapply(mycl, 1:nproj,
#                       function(obj) fasttab(D[,projs[,obj], drop=FALSE]))
    tabs <- parallel::parApply(mycl, projs,2,
                                function(obj) fasttab(D[, obj, drop=FALSE]))
    parallel::stopCluster(mycl); rm(mycl)
  }
  if (is.matrix(tabs)){
    ncovereds <- apply(tabs, 2, function(obj) sum(obj>0))
    whichnotcovereds <- apply(tabs, 2, function(obj) which(obj==0))
  }else{
    ## tabs is a list, because D has mixed levels
    ncovereds <- sapply(tabs, function(obj) sum(obj>0))
    whichnotcovereds <- lapply(tabs, function(obj) which(obj==0))
  }

    #### das war deutlich langsamer
     ### evtl. war diese Variante ohne parallel einen Hauch schneller
      # tots <- ncovereds <- rep(NA, nproj)
      # tabs <- whichnotcovereds <-
      #   vector(mode="list", length=nproj)
      #   doParallel::registerDoParallel(mycl)
      # foreach(i=1:nproj)  %do%{
      #   tots[i] <- prod(lls[projs[,i]])
      #   hilf <- fasttab(D[,projs[,i]])
      #   tabs[[i]] <- hilf
      #   ncovereds[i] <- sum(hilf>0)
      #   whichnotcovereds[[i]] <-
      #     which(hilf==0, arr.ind=TRUE)
      # }
  total <- sum(ncovereds)/sum(tots)
  proportions <- ncovereds/tots
  ave <- mean(proportions)
  min <- min(proportions)
  simple <- mean(ncovereds==tots)
  if (verbose==1)
    aus <- list(total=total, ave=ave, min=min, simple=simple,
                tots=tots, ncovereds=ncovereds, proportions=proportions)
  if (verbose==2)
    aus <- list(total=total, ave=ave, min=min, simple=simple,
                tots=tots, ncovereds=ncovereds,
                proportions=proportions,
                projections=projs,
                tabs=tabs,
                notcovereds=whichnotcovereds)
  ## remaining cases (anything but 1 and 2)
  if (!verbose %in% c(1,2)) aus <- list(total=total, ave=ave, min=min, simple=simple)
  class(aus) <- "coverage"
  aus
}

#' @importFrom gmp as.bigz
#' @importFrom gmp div.bigz
#' @export
coverage_iter <- function(D, t, isInteger=TRUE,
                     start0=TRUE, progress=FALSE, abortnot1=FALSE, start.proj=1){
  ## progress=TRUE costs a little bit of time for each projection
  stopifnot(is.matrix(D) || is.data.frame(D))
  if (is.data.frame(D) && isInteger) D <- as.matrix(D)
  if (is.matrix(D) && start0) D <- D + 1

  m <- ncol(D)
  if (t > m) stop("t > m")
  if (is.data.frame(D)){
    for (i in 1:m){
      D[[i]] <- as.integer(factor(D[[i]]))
    }
    D <- as.matrix(D)
  }
  ## now D is a matrix with values starting at 1
  ## (provided that the integer versions are adequate)
  lls <- levels.no.NA(D)
  if (start.proj==1 || !abortnot1){
      curproj <- 1:t
      start.proj <- 1
    }else{
    ## start at a later projection
    if (start.proj <= choose(m-1, t-1)){
      ## still in the projections that start with 1
      if (start.proj <= choose(m-2, t-2)){
          curproj <- 1:t
          start.proj <- 1
        }else{
          bounds <- cumsum(choose((m-2):(t-2),t-2))
          pos <- max(which(bounds <= start.proj))
          ## curproj is 1,(pos:(pos+t-2))+1
          curproj <- c(1, 2 + pos:(pos + t - 2))
          start.proj <- bounds[pos] + 1
      }
    }else{
      ## does not start with 1 any more
      bounds <- cumsum(choose((m-1):(t-1),t-1))
      pos <- max(which(bounds <= start.proj))
      curproj <- 1 + pos:(pos + t - 1)
      start.proj <- bounds[pos] + 1
    }
  }## now curproj is defined and can be iteratively increased

  nproj <- choose(m, t)

  ## initialize output objects

  simple <- gmp::as.bigz(0)  ## counts up and divides by nproj
  if (abortnot1 & start.proj>1) simple <- simple + start.proj - 1
  tot <- gmp::as.bigz(0)     ## counts up
  totcov <- gmp::as.bigz(0)  ## counts up and divides by tot
  ave <- 0     ## sums and divides by nproj
  if (abortnot1 & start.proj>1) ave <- ave + start.proj - 1
  min <- 1

  if (progress)
    pb <- txtProgressBar(0, nproj, initial=start.proj, char=":")

  tryCatch({
  for (i in start.proj:nproj){
    curtot <- prod(lls[curproj])
    curtab <- fasttab(D[,curproj,drop=FALSE])>0
    hilf <- sum(curtab)
    if (abortnot1) if (hilf < curtot){
      message("Coverage violated for projection ", paste(curproj, collapse=":"))
      aus <- list(total=NA,
                  ave=NA, min=NA,
                  simple=NA,
                  message=paste0("First projection with strength ", t,
                                 " coverage violated: ",
                                 paste(curproj, collapse=":"))
                  )
      class(aus) <- c("coverage", "list")
      return(aus)
    }
    simple <- simple + (hilf==curtot)
    tot <- tot + curtot
    totcov <- totcov + hilf
    ave <- ave + hilf/curtot
    min <- min(hilf/curtot, min)
    if (progress) setTxtProgressBar(pb, i)
    if (i<nproj) curproj <- gen.next.cbn(curproj, m)
  }
  if (progress) close(pb)
  aus <- list(total=as.numeric(gmp::div.bigz(totcov,tot)),
              ave=ave/nproj, min=min,
              simple=as.numeric(gmp::div.bigz(simple,nproj)))
  class(aus) <- c("coverage", "list")
  return(aus)  },
  interrupt = function(e) {
    message("Interrupted by user. Printing last completed projection id.")
    print(paste("Coverage ", t, " OK for first ", i-1, " projections."))
    if (abortnot1)
       print(paste("You might restart with start.proj set to ", i))
  })
}

prepcoverplot <- function(ts, type="projections",
                          main=NULL,
                          xlab=NULL, ylab=NULL,
                          las=1){
  ## should not be exported
  if (is.null(main)) main <- paste0("Coverage plot for t = ", paste(ts, collapse=", "))
  if (is.null(xlab)) xlab <- paste("%", type)
  if (is.null(ylab)) ylab <- "% covered"
  prepplot::prepplot(c(0,1), c(0,1), main=main, xlab=xlab, ylab=ylab,
           las=las, xaxs="i", yaxs="i", gridx=9, gridy=9,
           lwd.axis=1)

}

#' @export
coverplot <- function(D, t, isInteger=TRUE, start0=TRUE,
                      type="projections",
                      main=NULL,
                      xlab=NULL,
                      ylab=NULL, col=rgb(0,0,1,0.5), las=1,
                      plot=TRUE, ...){
  ## isInteger: are the values all from 0 to v-1 or from 1 to v?
  ## start0: 0 to v-1? (relevant for isI==TRUE only)
  hilf <- coverage(D, t, verbose=1, isInteger=isInteger, start0=start0)
  if (type=="projections"){
    tab <- table(hilf$proportions)
    x <- pctcombis <- c(0, rep(cumsum(rev(tab))/sum(tab), each=2), 0,0)
    y <- pctcovered <- c(rep(as.numeric(rev(names(tab))), each=2), 0,0,1)
  } else{
    tab <- table(hilf$proportions, hilf$tots)
    y <- rev(as.numeric(rownames(tab)))
    x <- rev(tab%*%as.numeric(colnames(tab)))
    x <- cumsum(x)/sum(x)
    x <- c(0, rep(x,each=2), 0,0)
    y <- c(rep(y,each=2),0, 0,1)
  }
  if(plot){
  prepcoverplot(t, main=main, xlab=xlab, ylab=ylab,
           type=type)
  polygon(x, y, col=col, ...)
  }
  invisible(list(x=x,y=y,t=t, coverage=hilf))
}
