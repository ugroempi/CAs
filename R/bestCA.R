#' Function to obtain the best uniform CA for a combination of t, k and v
#'
#' The function uses function Ns to find the best currently implemented
#' CA for the requested scenario.
#'
#' @rdname bestCA
#'
#' @importFrom utils download.file
#' @importFrom utils packageVersion
#' @importFrom curl has_internet
#'
#' @aliases bestCA
#' @aliases bestN
#' @aliases ckrsCA
#' @aliases dwyerCA
#' @aliases nistCA
#' @aliases tjCA
#'
#' @usage bestCA(t, k, v, fixNA=TRUE, seed=NULL, preference=NULL, override=FALSE, ...)
#' @usage bestN(t, k, v, internet=TRUE, exclude=NULL, ...)
#' @usage ckrsCA(t, k, v, ...)
#' @usage dwyerCA(t, k, v, ...)
#' @usage nistCA(t, k, v, ...)
#' @usage tjCA(t, k, v, ...)
#'
#' @param t integer: the strength
#' @param k integer: the number of columns
#' @param v integer: the number of levels for each column
#' @param fixNA logical: should flexible values be fixed?\cr
#'       If the ingoing CA has flexible values, \code{fixNA=TRUE}
#'       randomly assigns a fixed value to each flexible value.\cr
#'       Thus, a seed is needed for reproducibility.
#' @param seed \code{NULL} or an integer seed for making the
#'       fixing of NA values reproducible; if \code{seed=NULL},
#'       a seed is randomly obtained and stored in the attribute
#'       \code{fixNA_seed} of the returned object.
#' @param preference character: the label for the preferred method;
#'            one of "PALEY", "CAEX", "CYCLOTOMY", "CKRS",
#'            "recBoseCA_PCA", "recBoseCA_CA", "projBoseCA", "fuseBose",
#'            "WKS", "CS_MS", "CS_LCDST", "CS_CK", "DWYER", "NIST", "TJ"
#' @param override logical: If TRUE, the preference overrides a
#'           better method, otherwise (default) it only decides
#'           between equally-good methods.
#' @param internet logical; if FALSE, methods whose designs need an internet
#'           connection are excluded
#' @param exclude \code{NULL} or a character vector of method(s) to exclude,
#'   in support of size calculations for recursive constructions
#'   which need to calculate sizes of best ingredients
#' @param ... handed to construction function, currently,
#'          only implemented for \code{fuseBose}, for enabling \code{fixNA=FALSE}
#' @returns \code{bestN} returns a named integer: the first method that yields the
#'        smallest implemented size.\cr
#'        The other functions return a matrix of class \code{ca},
#'        which is a uniform covering array of strength at least \code{t} with
#'        \code{k} columns at \code{v} levels.
#'
#' @section Details:
#' Function \code{bestCA} uses function \code{\link{Ns}} for locating
#' the best construction for the given scenario, and subsequently
#' implements the first such construction, where the order is that
#' returned by function \code{\link{Ns}}.
#'
#' The user can override that order by stating a preferred construction
#' with the argument \code{preference}. With \code{override=TRUE},
#' the preference is observed whenever possible, even if this implies
#' a larger CA than necessary.
#'
#' Some constructions (\code{"DWYER"} and \code{"NIST"}) require an internet
#' connection, because the respective CAs are not part of the package.
#'
#' Some constructions provide CAs with flexible values, which are denoted
#' as NA values in package \pkg{CAs}. These need to be fixed for actual
#' experimentation, and function \code{bestCA} fixes them randomly; in such
#' cases, the returned object has an attribute \code{fixNA_seed}, which is
#' either the seed that was specified by the user or a seed that was created
#' by the function. Specifying the seed in a subsequent call to the function
#' permits to get the same CA again (provided that no new smaller CAs have
#' been implemented).
#'
#' The attributes of the returned object contain elements \code{eCAN},
#' \code{CAs-version} and \code{date}: \code{eCAN} gives the current
#' best array; as this may change, the other attributes provide the
#' context.
#'
#' @seealso [Ns()]
#' @examples
#' # example with a single best construction
#' Ns(2,12,4)
#' bestCA(2,12,4)
#' # example with three best constructions (TJ is now PALEY)
#' Ns(4,8,2)
#' # without a preference, bestCA uses Paley
#' bestCA(4,8,2)
#' bestCA(4,8,2, preference="CS_CK")
#'

#' Ns(3, 1500, 2)
#' ## powerCT is best
#' attributes(bestCA(3, 1500, 2))
#' ## it is, of course, also possible to directly use the
#' ## construction functions, in this case powerCA with type="CT".
#'
#' # a case that is not implemented
#' Ns(6, 45, 10)
#' try(bestCA(6, 45, 10))
#'
#' Ns(3, 50, 2)
#' ## DWYER needs internet connection, prefer TJ
#' D <- bestCA(3, 50, 2, preference="TJ")
#' dim(D)
#' # coverage(D, 3)
#'
#' Ns(6, 50, 2)
#' ## WKS is the only best
#' D <- bestCA(6, 50, 2)
#' dim(D)
#'
#' Ns(3,9,7)
#' ## CKRS and DWYER are best
#' D <- bestCA(3, 9, 7)
#' dim(D)
#'
#' Ns(3,9,6)
#' ## CKRS is best implemented, but 16 runs worse than eCAN
#' eCAN(3,9,6)
#' ## Torres-Jimenez designs are currently (June 2025) unavailable
#' D <- bestCA(3, 9, 6)
#' attributes(D)
#'
#' Ns(3,9,5)
#' ## CK_NRB is the unique best
#' ks(3,185,5)
#'
#' Ns(2,26,24) ## fuseBose is best-known
#' dim(bestCA(2,26,24))
#'
#' Ns(2,26,23) ## projBoseCA is best implemented with 623 runs
#' dim(D <- bestCA(2,26,23))
#' attributes(D)
#' ## eCAN 7 runs better from NCK post processing
#'
#' Ns(3,6,12)
#' ## the optimum Ji and Yin OA is in the package as oa1728.12.6
#'

#' @export
bestCA <- function(t, k, v, fixNA=TRUE, seed=NULL,
                   preference=NULL, override=FALSE, ...){
  Call <- sys.call()
  hilf <- Ns(t,k,v, ...)
  if (length(hilf)==1)
    stop("no construction for this setting has been implemented yet")
  ## return a matrix with v constant rows in case k <= t
  ## there will be an error otherwise
  if (k <= t) return(ffCA(k , v))
  internet <- curl::has_internet()
  if (!internet){
     message("there is no internet connection")
     if (!is.null(preference)) if (preference %in% c("DWYER","NIST"))
       stop(paste(preference, " requires an internet connection"))
     forbidden <- which(names(hilf) %in% c("DWYER","NIST"))
     if (length(forbidden) > 0)
        hilf <- hilf[-forbidden]
  }
  ## remove eCAN entry
  ## if it is feasible, there is another entry of the same size
  minN <- min(hilf[-length(hilf)])
  if (length(hilf)==0) stop("no feasible construction found")
  constrs <- names(hilf)[which(hilf<=minN)]
  ## handling of a stated preference
  finished <- FALSE
  if (!is.null(preference)){
    if (preference %in% constrs){
      aus <- eval(parse(text=labelToCode(preference, t, k, v, ...)))
      finished <- TRUE
    }
    else{
    message(preference, " is not among the best constructions.")
    if (override & preference %in% names(hilf)){
      message("It is used anyway because of override=TRUE.")
      aus <- eval(parse(text=labelToCode(preference, t, k, v, ...)))
      finished <- TRUE
    }else{
    if (preference %in% names(hilf))
      message("It is not used, because override=FALSE.")
    else
      message("It is not even among the possible constructions.")
    }}}
  if (!finished){
      aus <- eval(parse(text=labelToCode(constrs[1], t, k, v, ...)))
  }
  if (any(is.na(aus))){
    if (is.null(seed)) seed <- sample(1:32000, 1)
    nNA <- sum(is.na(aus))
    set.seed(seed)
    aus[is.na(aus)] <- sample(0:(v-1), nNA, replace=TRUE)
    attr(aus, "fixNA_seed") <- seed
  }
  dimnames(aus) <- NULL
  if (is.null(attr(aus, "origin"))) attr(aus, "origin") <- ifelse(finished, preference, constrs[1])
  if (is.null(attr(aus, "t"))) attr(aus, "t") <- t
  ## prepend the Call attribute
  attr(aus, "Call") <- c(Call, attr(aus, "Call"))
  attr(aus, "eCAN") <- eCAN(t,k,v)
  attr(aus, "date") <- Sys.Date()
  attr(aus, "CAs-version") <- packageVersion("CAs")
  if (!"ca" %in% class(aus)) class(aus) <- c("ca", class(aus))
  aus
}

#' @export
bestN <- function(t,k,v, internet=TRUE, exclude=NULL, ...){
  ## yields the size of the smallest implemented design,
  ## with or without internet
  stopifnot(is.numeric(t), is.numeric(k), is.numeric(v))
  stopifnot(t%%1==0, k%%1==0, v%%1==0)
  if (k <= t) return(c(FullFactorial=v^k))
  hilf <- Ns(t,k,v, exclude=exclude)
  if (internet) hilf <- hilf[setdiff(names(hilf), "eCAN")] else
    hilf <- hilf[setdiff(names(hilf), c("eCAN", "DWYER", "NIST"))]
  if (length(hilf)==0) return(NA)
  aus <- min(hilf)
  names(aus) <- names(hilf)[which.min(hilf)]
  aus
}

labelToCode <- function(label, t, k, v, ...){
  ## label must correspond to the label used in function Ns
  stopifnot(label %in% c("KSK","PALEY",
  "CAEX", "CYCLOTOMY", "CKRS", "miscCA", "recBoseCA_PCA", "SCA_Busht",
  "fuseBoseCA", "fuseBushtCA",
  "recBoseCA_CA", "projBoseCA", "compositCA",
  "WKS", "CS_MS",
  "CS_LCDST", "CS_CK", "powerCT", "DWYER", "NIST", "TJ", "CK_doublingCA",
  "CK_NRB","FullFactorial", "CS_CMMSSY"))
  if (label =="FullFactorial"){
    return(paste0("as.matrix(expand.grid(rep(list(0:(", v, "-1)),", k,")))"))
  }
  if (label =="KSK"){
    if (!v==2) stop('"KSK" requires v=2')
    return(paste0("KSK(", k, ")"))
  }
  if (label=="miscCA"){
    return(paste0("miscCA(", t, ", ", k, ", ", v, ", ...)"))
  }
  if (label=="compositCA"){
    return(paste0("compositCA(", t, ", ", k, ", ", v, ", ...)"))
  }
  if (label =="PALEY"){
    if (!v==2) stop('"PALEY" requires v=2')
    return(paste0("paleyCA(", t, ", ", k, ")"))
  }
  if (label =="fuseBoseCA"){
    return(paste0("fuseBoseCA(",k, ", ", v, ", ...)"))
  }
  if (label =="SCA_Busht"){
    return(paste0("SCA_Busht(",v,",", t, ")"))
  }
  if (label =="fuseBushtCA"){
    return(paste0("fuseBushtCA(",t, ",", k, ", ", v, ", ...)"))
  }
  if (label =="CAEX"){
    if (!v==3) stop('"CAEX" requires v=3')
    if (!t==2) stop('"CAEX" requires t=2')
    return(paste0("CAEX(", k, ")"))
  }
  if (label =="TJ"){
    ## as of June 2025, the 2-level CAs of TJ
    if (!v==2) stop('"TJ" requires v=2')
    return(paste0("tjCA(", t, ", ", k, ", ", v, ")"))
  }
  if (label =="CK_doublingCA"){
    if (!t==3) stop('"CK_doublingCA" requires t=3')
    return(paste0("CK_doublingCA(", k, ", ", v, ")"))
  }
  if (label =="CK_NRB"){
    if (!t==3) stop('"CK_NRB" requires t=3')
    if (!v %in% 3:5) stop('"CK_NRB" requires v = 3, 4, or 5')
    if (!k<=2*v) stop('"CK_NRB" requires k<=2*v')
    return(paste0("CK_NRB(", k, ", ", v, ")"))
  }
  if (label=="CYCLOTOMY"){
    return(paste0("cyclotomyCA(", t, ", ", k, ", ", v, ")"))
  }
  if (label=="CKRS"){
    return(paste0("ckrsCA(", t, ", ", k, ", ", v, ")"))
  }
  if (label=="recBoseCA_PCA"){
    if (!t==2) stop('"recBoseCA_PCA" requires t=2')
    if (!v %in% primedat$q) stop("v must be prime or prime power")
    return(paste0("recBoseCA(", t, ", ", k, ", ", v, ", type='PCA')"))
  }
  if (label=="recBoseCA_CA"){
    if (!t==2) stop('"recBoseCA_CA" requires t=2')
    if (!v %in% primedat$q) stop("v must be prime or prime power")
    return(paste0("recBoseCA(", t, ", ", k, ", ", v, ", type='CA')"))
  }
  if (label=="projBoseCA"){
    if (!t==2) stop('"projBoseCA" requires t=2')
    if (!v %in% primedat$q) stop("v must be prime or prime power")
    return(paste0("projBoseCA(", k, ", ", v, ")"))
  }
  if (label=="WKS"){
    if (!t==6) stop('"WKS" requires t=6')
    if (!v==2) stop('"WKS" requires v=2')
    ## the CAs are labeled by k, and there is an entry
    ## for each k between 30 and 72 (including)
    return(paste0('WKS_CAs[["', k, '"]]'))
  }
  if (label=="CS_MS"){
    if (!t==2) stop('"CS_MS" requires t=2')
    return(paste0('CS_MS(', k, ', ', v, ')'))
  }
  if (label=="CS_LDCST"){
    if (!t==2) stop('"CS_LCDST" requires t=2')
    return(paste0('CS_LCDST(', k, ', ', v, ')'))
  }
  if (label=="CS_CMMSSY"){
    if (!t==2) stop('"CS_CMMSSY" requires t=2')
    return(paste0('CS_CMMSSY(', k, ', ', v, ')'))
  }
  if (label=="CS_CK"){
    if (!v==2) stop('"CS_CK" requires v=2')
    return(paste0('CS_CK(', k, ', t=', t, ')'))
  }
  if (label=="powerCT"){
    return(paste0('powerCA(', t, ', ', k, ', ', v, ', type="CT")'))
  }
  if (label=="DWYER"){
    return(paste0('dwyerCA(', t, ', ', k, ', ', v, ')'))
  }
  if (label=="NIST"){
    return(paste0('nistCA(', t, ', ', k, ', ', v, ')'))
  }
}

#' @export
ckrsCA <- function(t, k, v, ...){
  Call <- sys.call()
  hilf <- CKRScat[CKRScat$t>=t & CKRScat$v==v & CKRScat$k>=k,]
  if (nrow(hilf)==0) stop("there is no suitable array in CKRS_CAs")
  N <- min(hilf$N)
  tachieved <- hilf$t
  hilf <- hilf[which.min(hilf$N),,drop=FALSE]
  N <- hilf$N
  fn <- hilf$fn
  aus <- CKRS_CAs[[fn]]
  attr(aus, "t") <- tachieved
  attr(aus, "Call") <- Call
  aus[,1:k]
}

#' @export
dwyerCA <- function(t, k, v, ...){
  Call <- sys.call()
  hilf <- DWYERcat[DWYERcat$t>=t & DWYERcat$v==v & DWYERcat$k>=k,]
  if (nrow(hilf)==0) stop("there is no suitable array in the Dwyer database")
  N <- min(hilf$N)
  hilf <- hilf[which.min(hilf$N),,drop=FALSE]
  N <- hilf$N  ## should not be necessary, just to be safe
  kneeded <- hilf$k
  tachieved <- hilf$t
  nam <- paste0(paste("CA", N, t, kneeded, v, sep="_"), ".txt")
  ## load from Dwyer repo
  pfadGithub <- "https://raw.githubusercontent.com/aadwyer/CA_Database/main/Repository/CA"
  pfad <- paste0(pfadGithub, "/", nam)
  aus <- readCA(pfad, ninstruct = 0, ignore.chars=c("[","]"),
                header=FALSE, sep=",",
                origin=paste0("Dwyer Github repository, ", nam, ", ",
                              hilf$Source))
  attrs <- attributes(aus)
  attrs$dimnames <- NULL
  aus <- aus[,1:k]
  attrs$dim <- dim(aus)
  N <- nrow(aus)
  neles <- lengths(lapply(1:N, function(obj) unique(aus[obj,])))
  constrows <- which(neles==1)
  if (length(constrows) > 1) aus <- aus[c(constrows, setdiff(1:N, constrows)),]
  attributes(aus) <- attrs
  attr(aus, "t") <- tachieved
  if (length(constrows) > 1) attr(aus, "nconstant") <- length(constrows)
  attr(aus, "Call") <- Call
  aus
}

#' @export
nistCA <- function(t, k, v, ...){
  Call <- sys.call()
  hilf <- NISTcat[NISTcat[,"t"]==t & NISTcat[,"v"]==v &
                    NISTcat[,"k"]==k,, drop=FALSE]
  if (nrow(hilf)==0) stop("there is no suitable array in the NIST library")
  nam <- rownames(hilf)
  pfad <- paste0("https://math.nist.gov/coveringarrays/ipof/cas/t=", t,"/v=",v,"/", nam)
  speicherpfad <- paste0(tempdir(), "/", nam)
  aus <- try(download.file(url = pfad, destfile = speicherpfad))
  if ("try-error" %in% class(aus)) message("error encountered in downloading the array")
  if (!"try-error" %in% class(aus)){
    txtpfad <- gsub(".zip","",speicherpfad, fixed=TRUE)
    utils::unzip(speicherpfad, exdir = tempdir())
    aus <- readCA(txtpfad, ninstruct=0, skiplines = 1)
  }
  attr(aus, "origin") <- "NIST covering array library"
  attr(aus, "t") <- t
  attr(aus, "Call") <- Call
  aus
}

#'@export
tjCA <- function(t, k, v=2, ...){
  Call <- sys.call()
  if (!v==2) stop("For v=3, use CAEX. v>3 was not yet implemented.")
  hilf <- TJcat[TJcat$t>=t & TJcat$v==v & TJcat$k>=k,]
  if (nrow(hilf)==0) stop("there is no suitable array in TJ2level_CAs")
  N <- min(hilf$N)
  hilf <- hilf[which.min(hilf$N),,drop=FALSE]
  N <- hilf$N  ## should not be necessary, just to be safe
  kneeded <- hilf$k
  tachieved <- hilf$t
  nam <- hilf$nameInTJ2level_CAs
  if (!nam==""){
  aus <- TJ2level_CAs[[nam]][,1:k]
  attr(aus, "origin") <- "Torres-Jimenez repository, Feb 6 2025"
  attr(aus, "t") <- tachieved
  attr(aus, "Call") <- Call
  return(aus)
  }
  code <- hilf$code
  if (code=="") stop("this setting is not implemented in tjCA, \nuse function Ns for finding alternatives")
  aus <- eval(parse(text=code))
  attr(aus, "comment") <- "called via tjCA"
  aus
}

ffCA <- function(k, v){
  Call <- sys.call()
  aus <- maxconstant(as.matrix(expand.grid(rep(list(0:(v-1)), k))))
  attr(aus, "t") <- NA
  attr(aus, "origin") <- "Full Factorial"
  attr(aus, "nconstant") <- v
  attr(aus, "Call") <- Call
  aus
}
