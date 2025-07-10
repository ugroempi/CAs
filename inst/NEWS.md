# NEWS

10 July 2025, version 0.8
 - implemented data file powerCTcat in support of power
   construction according to Colbourn and Torres-Jimenez (2010)
 - implemented functions powerCA (with type="CT"), N_powerCT and 
   k_powerCT and incorporated them into Ns, bestN, bestCA
 - removed seven arrays from TJcat and TJ2level_CAs, as well 
   as from the folder TJ2level in extdata
 - moved SCA_Bose and SCA_Busht to a common file with dedicated
   documentation
 - fixed various bugs in other constructions that showed up in 
   testing the new function powerCA

05 July 2025, version 0.7
 - implemented a data file miscCAcat, a function miscCA 
   with N_ and k_ functions,
   and modified Ns and ks to refer to these
 - modified oa1728.12.6 to PCA structure, and added further 
   cas / oas, as well as retrieval of a few OAs from DoE.base
   (handled via miscCAcat);
   (renamed the R code in extdata for single arrays to miscCAs.R)
 - implemented the three arrays of source "Chateauneuf-Kreher NRB"
 - rewrote and re-exported k_fuseBose
 - removed TJcat arrays that can be replaced by other constructions, 
   and adapted TJcat so that they can still all be constructed
 - bug fix to productPCA: did not work as expected when trying 
   without success to improve a PCA with tryhard=TRUE

02 July 2025, version 0.6

 - added a function bestN for use in size calculations
   for recursive constructions
 - added an internet availability check to function bestCA
 - modified ckrsCA, nistCA, dwyerCA and tjCA to throw 
   a meaningful error for cases that are not covered
 - extended CS_LCDST to cases with f>3 
   (uses bestCA for a suitable ingredient)
 - fixed LCDSTStarters for v=10 and k=27 and
   for v=11 and k=26, and correspondingly
   modified LCDSTCombis
 - modified LCDSTCombis$N to hold realistic
   numbers for f>3 based on the currently available
   minimum size construction (three different columns for N now)
 - added a forgotten entry to DWYERcat, and added 
   useful row names 
 - added folder extdata with all data files that 
   are used by the package and not available 
   from the internet, as well as code for creating 
   the objects in data folder
 - deleted function CK_constrD.R, 
   which was not really for construction D but a duplicate
   of CK_doubling
 - added an argument "exclude" to functions Ns and bestN, 
   in support of recursive constructions that need best sizes
   for ingredients
 - immediately return Bose CA for small k with CAEX, 
   as lineage for N=9 does not exist
 - automated CK_doubling via CK_doublingCA with arguments k and v,
   added an N-function for it and included it in Ns and bestCA
 - bug fix for documentation of DHHF

30 June 2025, version 0.5

 - bug fix: prevented nonsensical result from N_BoseCA 
 - bug fix: fixed mistake in N_fuseBose and N_fuseBose_forNs
 - bug fix: fixed mistake in LCDSTCombis (for cases with NAs, 
   f was one too large, resulting in too small N; occurred for 7 rows;
   led to erroneous N reported, but arrays were correct)
 - bug fix: fixed mistake in LCDSTStarters for v=7 with k=14
 - bug fix: CS_LCDST erroneously refused creation of an array in case of 
   f=3 with NA-values in the starter
 - added function fuseBoseCA, in line with the other constructions
 - renamed N_fuseBose_forNs to N_fuseBoseCA, 
   in line with the other constructions
 - bug fix: hid k_fuseBose, because fixing it takes more care (and it is 
   not so urgently needed)

27 June 2025, version 0.4

 - the user-visible data are now in a data folder instead of sysdata.rda, 
   and are documented in more detail
 - the Colbourn and Keri (2009) construction was implemented in function CS_CK,
   and also added to the repertoire for bestCA
 - added fuseBose and SCA_Bose to Ns and bestCA
 - increased the version number

24 June 2025, Version 0.3
 
 - implemented function bestCA for obtaining the best CA from 
   the currently implemented constructions, based on function Ns
 - implemented functions ckrsCA, dwyerCA, nistCA 
   and tjCA for creating designs from catalogues
 - added an internal utility function labelToCode that 
   created the suitable code line for function bestCA
 - fixed a bug in Ns that led to a deep recursion
 - adapted the print method for class ca to new attributes
 - added an argument skiplines to readCA
 

22 June 2025, still Version 0.2.2

 - renamed PHF2CA to DHF2CA, as a PHF is a special case of a DHF
 - added TJ2level_CAs in support of CT constructions
 - made TJcat into data.frame with info column on element name for TJ2level_CAs
 - implemented DHHF2CA for heterogeneous DHF CT construction
 - added oa1728.12.6 (OA of strength 3 with six 12-level columns based on 
   Ji/Ying construction); this may eventually be moved supplemented with 
   other Ji/Ying OAs, and moved to DoE.base

18 June 2025, still Version 0.2.2

 - renamed file CK.R to CK_constrD.R, in line with its content and its documentation file
 - fixed a bug in projectionBose (two many rows were made constant)
 - fixed a bug in cyclotomyCA (it always returned the maximum number of columns instead of 
   the requested k columns)
 - improved PHF2CA to better take advantage of constant rows, in the sense of corollary 2.4
   of Colbourn and Torres-Jimenez 2010

05 June 2025, still Version 0.2.2
 
 - renamed previous function projectionBose to fuseBose, as it does a specific type
   of fusing in line with Colbourn et al. (2010) rather than a projection
 - added new function projectionBose that does indeed obtain a projection
   from a Bose array (Theorem 2.3 of Colbourn 2008; output CA can be uniform or mixed level)
 - added new function projBoseCA that allows to input k and v and obtain a projection 
   based uniform CA (uses projectionBose as workhorse function)

25 May 2025, still Version 0.2.2
 
 - added function kd for calculating the k for recursiveBose constructions
   of type PCA (in some instances based on experimental results), 
   exported this function and function Dd
 - made function recursiveBose (and thus recBoseCA) fully exploit the 2PCA 
   structure of the SCA_Bose array
 - did the same for function productPCA, if D2=NULL

22 May 2025, Version 0.2.2

 - added a function PHF2CA, in support of later implementation of Power 
   constructions
 - added a function recBoseCA for creating CAs from recursive processing 
   of several SCA_Bose CAs; 
   that recursive processing was changed from the method by Hartman 2005 
   (which is still there as recursiveBoseHartman) to 
   using either productPCA (which has the same number of rows but more 
   columns as the Hartman method) 
   or productCA (which has both more rows and more columns), 
   governed by a type argument in functions recBoseCA and recursiveBose
 - added further Ns and ks functions for post-processing methods
 - changed readCA to take a single path argument including the file name
 - added names to the entries of list CKRS_CAs
 - bug fix for documentation of productCA and productPCA: 
   run size for D2 was equal to run size for D1, changed it from N to M
 - BUG FIX: productCA with generalized=TRUE sometimes failed to make the 
   constant rows distinct, which led to violation of coverage (CAEX was not affected).


19 May 2025, still Version 0.2.1

 - added functions Ns and ks
 - modified N and k functions, so that most have compatible arguments
   renamed prior function k_CAEX to k_detail_CAEX
 - added names to list entries in WKS_CAs
 - added data.frame CYCLOTOMYcat and functions cyclotomyCA, N_CYCLOTOMYcat and k_CYCLOTOMYcat
 - augmented the internal object primedat with a few more primes in support of CYCLOTOMYcat


17 May 2025, still Version 0.2.1

 - added Paley constructions according to Colbourn (2015)
 

14 May 2025, still Version 0.2.1

 - made function `CK_doubling` use the `CAEX` designs per default for 3-level cases
   (i.e., `D2=NULL` now works for v=2 and v=3)
 - updated function `postopNCK` from a very naive version to a more serious implementation
   of the algorithm (which has still not been successful in any challenging case)
 - added functions `crossSum` and `directSum`
 - modified `readCA` to also read data without separators, where levels 10, 11, ... 
   are coded as lower case letters a, b, ...
 - modified `coverage` to also work with `t=1`
 - added data.frame `CKRScat` and list `CKRS_CAs`
 - added several query functions: `N_DWYERcat`, `N_WKScat`, `N_CKRScat`, 
   `k_DWYERcat`, `k_WKScat`, `k_CKRScat`


10 May 2025: Version 0.2.1

 - added function `CS_LCDST` for strength 2 cover starter constructions of Lobb et al. 2012,
   with internal list `LCDSTStarters` and internal data.frame object `LCDSTCombis` with 
   information on what these construct. 
 - exported function `SCA_MS` for a rearranged Meagher/Stevens Cover Starter array.
 - added function `projectionBose` to construct array in q^2-3 runs for q+1 factors
   in q-1 levels (q prime power), according to Colbourn 2008, with corresponding N- 
   and k-functions
 - implemented `postopNCK` to remove rows according to the proposal of Nayeri, Colbourn 
   and Konjevod (2013)
 - added internal data.frame `DWYERcat` with information about the DWYER (2024) CAs database
   (not yet in any productive use)
 - Bug fix: `Hartman74` returned fewer than `k^2` columns, due to unnecessary removal of 
   constant rows in `B`

05 May 2025: Version 0.2.0 

  - renamed function `productCA1` to `productCA` and function `productCA2` to `crossCAs`;
      separated the documentation file `productCA.Rd` into three files, corresponding 
      to function names
  - the direct product function `productCA` now also yields the *generalized* direct 
      product, in the sense that it omits at least 2 and up to v constant rows from 
      the ingoing CAs (after making rows constant, if necessary)
  - there are now also functions to calculate k as a consequence of N for the 
      libraries of designs (i.e., `eCAK`, `k_NISTcat`, `k_TJcat`)
  - `CAEX` gained a logical argument `maxk1`, and the internal object `CAEX_lineages`
      gained the information from Torres-Jimenez et al. (2021) about the achievable k1, 
      which is used as the target for k1 maximization.
  - various improvements to documentation
  - BUG fix: latest state of version 0.1.0 had a serious bug in `productPCA` 
      and consequently in the PCA-based instances of `CAEX`, 
      which caused violation of coverage in most situations.
      This was fixed by using the new internal function `productCA_raw` 
      within `productPCA`.

Version 0.1.0: first version on GitHub
