## Colbourn table rows for strength t and v levels, augmented with the best N
## (and its construction) from package CAs, plus extra rows for k values at
## which bestN changes but eCAN (Colbourn table) does not.
##
## Cheap route: the Colbourn table is a step function in k with breakpoints at
## its listed k values. bestN is evaluated only at these breakpoints, and
## within a stretch of constant eCAN, further change points of bestN are found
## by bisection. Relies on: bestN is non-decreasing in k, so that
## equal values at both ends of a stretch imply a constant stretch.
##
## Extra rows are labelled by the LAST k of each plateau of bestN, i.e., by the
## k with bestN(k) < bestN(k+1), analogous to the maximum k per N in the table.
## Initial code was created by Claude 5.5

library(CAs)
bestN_table <- function(t = 3, v = 3, kmin = 4, kmax = 10000, ...){
  cb <- CAs::colbournBigFrame
  cb <- cb[cb$t == t & cb$v == v & cb$k >= kmin & cb$k <= kmax, ]
  cb <- cb[order(cb$k), ]
  cb <- cb[!duplicated(cb$k), ]   # one row per k
  stopifnot(nrow(cb) > 0)
  cache <- new.env()
  f <- function(k){
    key <- as.character(k)
    hit <- cache[[key]]
    if (is.null(hit)){
      b <- bestN(t, k, v, ...)
      hit <- list(N = if (length(b) == 0 || is.na(b)) NA_real_ else unname(b),
                  constr = if (length(b) == 0 || is.na(b)) NA_character_
                  else names(b))
      assign(key, hit, envir = cache)
    }
    hit
  }
  Nval <- function(k){ n <- f(k)$N; if (is.na(n)) Inf else n }
  ## plateau ends k in [lo, hi-1] with bestN(k) < bestN(k+1),
  ## given bestN(lo) <= ... <= bestN(hi)
  changes <- function(lo, hi){
    if (lo >= hi || Nval(lo) == Nval(hi)) return(integer(0))
    if (hi - lo == 1) return(lo)
    mid <- (lo + hi) %/% 2
    c(changes(lo, mid), changes(mid, hi))
  }
  ## stretch i covers the k values prev+1, ..., cb$k[i], on which eCAN is cb$N[i]
  starts <- c(kmin, head(cb$k, -1) + 1)
  extra <- unlist(lapply(seq_len(nrow(cb)),
                         function(i) changes(starts[i], cb$k[i])))
  extra <- setdiff(extra, cb$k)
  ks <- sort(unique(c(cb$k, extra)))
  eC <- lapply(ks, function(k) eCAN(t, k, v))
  aus <- data.frame(
    k = ks,
    eCAN = (ref <- vapply(eC, function(obj) obj$CAN, numeric(1))),
    Source = vapply(eC, function(obj) as.character(obj$Source), character(1)),
    bestN = (implem <- vapply(ks, function(k) f(k)$N, numeric(1))),
    construction = vapply(ks, function(k) f(k)$constr, character(1)),
    diff = implem - ref,
    'row type' = ifelse(ks %in% cb$k, "Colbourn", "extra"))
  aus
}
## example (t=3, v=3, k from 4 to 10000); further arguments go to bestN
# tab <- bestN_table(3, 3, 4, 10000)
# head(tab, 20)
# table(tab$extra)

for (v in 2:25){
  for (t in 2:6){
    tab <- bestN_table(t, v, t+1, 10000, maxfuse=6)
    require(gt)
    gttab <- gt(tab)
    gttab <- tab_style(
      data = gttab,
      style = cell_fill(color = "lightgreen"), ## #90EE90;
      locations = cells_body(columns = diff, rows = diff <= 0)
    )
    gttab <- tab_style(
      data = gttab,
      style = cell_fill(color = "lightyellow"),  ## #FFFFE0;
      locations = cells_body(columns = diff, rows = diff > 0 & diff <= eCAN*0.05)
    )
    gttab <- tab_style(
      data = gttab,
      style = cell_fill(color = "orange"),  ## #FFA500;
      locations = cells_body(columns = diff, rows = diff > eCAN*0.1)
    )
    gttab <- tab_spanner(gttab, "best N from package CAs 0.25", 4:5)
    gttab <- tab_spanner(gttab, "eCAN from Colbourn tables", 2:3)
    gttab <- tab_header(gttab, paste0("Colbourn tables with CAs bestN values for t=", t, " and v=", v),
                        subtitle="diff color coded: green=bestN equal or better, yellow=bestN up to 5% worse, orange=bestN more than 10% worse")
    gtsave(data=gttab, filename = paste0("t",t,"v",v,".html"),
           path="c:/rtests/CAsbak/bestNtables/")
  }
}
