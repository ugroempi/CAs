#' coverage and diagnostic plot functions
#'
#' functions for calculating covering proportions for a given strength,
#' and for visualizing those
#'
#' @rdname coverage
#' @aliases coverage
#' @aliases prepcoverplot
#' @aliases coverplot
#'
#' @usage coverage(D, t, isInteger=TRUE, verbose=0, start0=TRUE, parallel=1)
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
#' with increasing detail added for verbose 1 or 2.
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
#' If verbosity is requested, there are additional list elements:
#'
#' for \code{verbose=1} the number of covered tuples for each projection (\code{ncovereds}),
#' the proportion this number represents (\code{proportions}),
#'
#' for \code{verbose=2} the frequency tables of tuples (\code{tabs})
#' and the identification indices of uncovered tuples \code{uncovereds}.
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
  ## start0: 0 to v-1? (relevant for isI==TRUE only)
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
  if (is.matrix(D) && start0) D <- D + 1

  m <- ncol(D)
  if (is.data.frame(D)){
    for (i in 1:m){
      D[[i]] <- as.integer(factor(D[[i]]))
    }
    D <- as.matrix(D)
  }
  ## now D is a matrix with values starting at 1
  ## (provided that the integer versions are adequate)
  lls <- levels.no.NA(D)
  projs <- nchoosek(m,t)
  nproj <- ncol(projs)
  ## make D into integer matrix from 1 to lls

  tots <- sapply(1:nproj, function(obj) prod(lls[projs[,obj]]))
  if (parallel==1)
    tabs <- lapply(1:nproj, function(obj) fasttab(D[,projs[,obj], drop=FALSE]))
  else{
    stopifnot(requireNamespace("parallel"))
    stopifnot(parallel <= parallel::detectCores())
    mycl <- parallel::makePSOCKcluster(parallel)
    tabs <- parallel::parLapply(mycl, 1:nproj,
                       function(obj) fasttab(D[,projs[,obj], drop=FALSE]))
    parallel::stopCluster(mycl); rm(mycl)
  }
  ncovereds <- sapply(tabs, function(obj) sum(obj>0))
  whichnotcovereds <- lapply(tabs, function(obj) which(obj==0))

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
