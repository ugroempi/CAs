#' Construction Chateauneuf-Kreher NRB of strength 3 for v=3,4,5 with k=2*v
#'
#' Function for three very good strength 3 CAs with v=3,4,5 and k=2*v
#' from Chateauneuf, Colbourn and Kreher 1999 and Chateauneuf and Kreher 2002
#'
#' @rdname CK_NRB
#'
#' @aliases CK_NRB
#' @aliases N_CK_NRB
#' @aliases k_CK_NRB
#'
#' @usage CK_NRB(k, v)
#' @usage N_CK_NRB(t=3, k, v)
#' @usage k_CK_NRB(t=3, N, v)
#'
#' @param t integer, strength (must be 3)
#' @param k integer, number of columns needed (<= 2*\code{v})
#' @param v integer, number of levels for each column
#' @param N integer, maximal number of affordable runs
#'
#' @returns \code{CK_NRB} returns a covering array of 33 (\code{v}=3),
#' 88 (\code{v}=4) or 185 rows (\code{v}=5),\cr
#' \code{N_CK_NRB} returns a number for valid entries and \code{NA} otherwise.
#'
#' @section Details:
#' CA(33,3,6,3) is from Figure 1 of Chateauneuf, Colbourn and Kreher (1999),\cr
#' CA(88,3,8,4) is from Example 2.3 of the paper with the ingredient M for
#' from Figure 2 (where typos in the last row of M had to be fixed).\cr
#' CA(185,3,10,5) is from Theorem 3.4 of Chateauneuf and Kreher (2002,
#' with AGL(2,q) as stated above Theorem 3.1 and the ingredient M in Figure 5 of the paper).
#'
#' @references Chateauneuf, Colbourn and Kreher (1999), Chateauneuf and Kreher (2002)
#'
#' @examples
#' dim(CK_NRB(6, 3))
#' attributes(CK_NRB(9,5))
#' Ns(3, 7, 4)  ## substantially better than other constructions
#' Ns(3, 9, 4) ## not possible
#' dim(bestCA(3,8,4))
#' k_CK_NRB(3, 100, 4) ## with up to 100 runs, 8 columns affordable
#'

#' @export
CK_NRB <- function(k, v){
  stopifnot(v %in% 3:5)
  stopifnot(k <= 2*v)
  ## the function is dealt with in three portions
  if (v==3){
      ## copy-paste from Chateauneuf Colbourn Kreher 1999
      aus <- t(CAs:::funmakefromstrings(
      c("012211200220110021121022021001012",
      "122102002101102211200220110012012" ,
      "221010021211020112022201000121012" ,
      "210120212010201120212010201210012" ,
      "101222120002011202110102212100012" ,
      "000001111122222000001111122222012")))[,1:k]
      class(aus) <- c("ca", class(aus))
      attr(aus, "t") <- 3
      attr(aus, "origin") <- "Chateauneuf Colbourn Kreher 1999, Figure 1"
      attr(aus, "Source") <- "Chateauneuf-Kreher NRB"
      return(aus)
  }
  if (v==4){
      ## creation of ca88.4.8
      ## M is copy-paste from Chateauneuf Colbourn Kreher 1999
      ## with typos fixed
      M <- CAs:::funmakefromstrings(
      c("0000000",
      "0233121" ,
      "1023312" ,
      "2102331" ,
      "1210233" ,
      "3121023" ,
      "3312102" ,
      "2331210"))  ## 8 x 7 typos fixed
      ## "2221210"))  ## 8 x 7 ## fixed

      S4 <- permutations::allperms(4)
      A4 <- S4[permutations::is.even(S4)]
      A4fun <- permutations::as.function.permutation(A4)
      A4funvec <- Vectorize(A4fun)
      ## each column is the function applied to one element of the matrix
      ## each row is a matrix, in column-wise order
      mat <- matrix(0:3, 8,4, byrow=TRUE)
      hilf <- A4funvec(M+1) - 1
      mat2 <- do.call(cbind, lapply(1:nrow(hilf), function(obj)
        matrix(hilf[obj,], nrow=8)))
      aus <- t(cbind(mat,mat2))[,1:k]
      class(aus) <- c("ca", class(aus))
      attr(aus, "t") <- 3
      attr(aus, "origin") <- c("Chateauneuf Colbourn Kreher 1999, Example 2.3",
                                    "based on corrected Figure 2")
      attr(aus, "Source") <- "Chateauneuf-Kreher NRB"
      return(aus)
  }
  if (v==5){
        ## ca185.5.10 of Theorem 3.4 of Chateauneuf & Kreher 2002
        ## who reference unpublished work of Kreher and Radzisowski for M
        ## the only M that does it, according to Chateauneuf and Kreher
        ## (typed in from Chateauneuf and Kreher 2002)
        M <- rbind(
          rep(0,9),
          c(0,rep(1,8)),
          c(1,1,0,2,4,2,3,3,4),
          c(1,0,1,4,2,3,2,4,3),
          c(2,4,3,1,0,2,4,2,3),
          c(2,3,4,0,1,4,2,3,2),
          c(4,3,2,4,3,1,0,2,4),
          c(4,2,3,3,4,0,1,4,2),
          c(3,4,2,3,2,4,3,1,0),
          c(3,2,4,2,3,3,4,0,1)
        )

        # table(M) ## sanity check

        ## AGL(2, q) is the set of functions f: x ->  (ax + b)%%q
        ##                              with a,b from GF(q) and a non-zero

        # list of the elements of AGL(2,5) consists of this function
        # applied with (a,b) from variants
        variants <- expand.grid(a=1L:4L, b=0L:4L)

        lapply(1:nrow(variants), function(obj)
               function(x) (variants[obj,]$a*x+variants[obj,]$b)%%5)

        ## initially in transposed format
        ## the constant rows
        mat <- matrix(0:4, 10, 5, byrow=TRUE)
        ## apply the 20 functions (including a=1 and b=0 for identity)
        mat2 <- do.call(cbind, lapply(1:nrow(variants),
                                      function(obj) matrix(AGL2(M, 5, variants[obj,]$a,
                                                            variants[obj,]$b), nrow=10)))
        aus <- t(cbind(mat, mat2))[,1:k]
        class(aus) <- c("ca", class(aus))
        attr(aus, "t") <- 3
        attr(aus, "origin") <- c("Chateauneuf Colbourn Kreher 2002, Theorem 3.4",
                                 "M typed in")
        attr(aus, "Source") <- "Chateauneuf-Kreher NRB"
        return(aus)
  }
}

#' @export
N_CK_NRB <- function(t=3, k, v){
  if (!t==3) return(NA)
  if (!v %in% c(3,4,5)) return(NA)
  if (!k<=2*v) return(NA)
  switch(as.character(v), NA, '3'=33, '4'=88, '5'=185)
}

#' @export
k_CK_NRB <- function(t=3, N, v){
  if (!t==3) return(NA)
  if (!v %in% c(3,4,5)) return(NA)
  if (!N>=switch(as.character(v), NA, '3'=33, '4'=88, '5'=185)) return(NA)
  switch(as.character(v), NA, '3'=6, '4'=8, '5'=10)
}
