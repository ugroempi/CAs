# NEWS

 - TODO: check Call attribute for MCA2 and postopNCK, in the combination
   it apparently loses the MCA2 call 
   (see file D15etalOptimized.RData in colbournTableReport folder)
 - TODO: check handling of flexible values by postopNCK and the mixed level 
   design creation functions
   
February 21 2025, version 0.23
 - added functionality for SCAs

November 21 2025, version 0.22
 - added arrays from Walker and Colbourn 2009 for scphfCA
 - added classical (poor) augmentation constructions 
   (CK_constrD and CKRS_augment), not included in function Ns
 - added the L81.3.5 to miscCAcat
   
October 11 2025, version 0.21
 - improved the implementation of fusion in bestN, so that the speed of bestN 
   is again in its prior range
 - added function add1
 - added a reference for crossCAs
   
September 27 2025, version 0.20
 - improved N_CS_LCDST to return the minimum and not the first entry, because 
   a later entry is occasionally better
 - improved MCA2 so that the cheap reduction of run size is always included, even 
   with outerRetry=0 (in auxiliary function outerseparate)
 - modified MCAt to set levels to missing rather than that of last row 
   (improves chances for immediate row removal for easy cases)
 - improved postopNCK to remove rows that have fewer than t non-flexible values
 - updated CPHFcat and WCS_CAs to updated source status, and updated powerCTcat 
   accordingly

September 23 2025, still version 0.19
 - improved N_powerCT by restricting the search to strengths t and t+1
 - added a maxN argument to bestCA

September 22 2025, version 0.19
 - automated CA_to_MCA in new function MCAt
 - implemented CPHFs of Wagner, Colbourn and Simos (2022)
 - reduced powerCTcat to constructions with at most 1000000 runs
 - implemented fuse for non-prime power v for bestN and bestCA, 
   and removed the preference and override arguments from bestCA in the process
 - fixed handling of flexible values in productCA and CAEX
 - improved postopNCK: it now uses a good optimality bound for t>2 as well as for t=2

September 05 2025, version 0.18
 - function pcaCA based on the data frame PCAcat:
   implemented productPCA constructions with promising PCA ingredients (large k1)
   from Bose (not combining only Bose, because these are in recursiveBose), all cover starter 
   constructions with only one fixed value, and some miscCA designs
 - function dpCA based on the data frame DPcat:
   implemented productCA constructions with promising CA ingredients (constant rows)
   from Bose (not combining only Bose, because these are in recursiveBose), all cover starter 
   constructions with only one fixed value, and some miscCA designs; (basically copied the list 
   from productPCA, more might be achievable with more deliberate array selection)
 - implemented productPCA constructions with promising PCA ingredients (large k1)
   from Bose (not combining only Bose, because these are in recursiveBose), all cover starter 
   constructions with only one fixed value, and some miscCA designs
 - adapted LCDSTCombis with some improved ingredients
 - adapted power constructions with improved ingredients - apparently nothing changed, 
   but this can not be verified because the old rda file on GitHub appeared to be corrupt
 - added arguments c1 and c2 for productCA, stopped it from enforcing a PCA (it now instead 
   keeps its constant rows at the top), and changed the default for check to FALSE
 - included the Cohen simulated annealing CAs and two CAs from Kokkala et al. (2018) into miscCA
 - included an array from the paper CMMSSY (2006) into miscCA
 - modified PCAstatus to k in miscCAcat for arrays with v constant rows, 
   and sorted miscCAcat by t, v and k
 - removed CA_to_PCA from productCA and added it in in CAEX
 - added new arguments for function productCA with the intention of making it faster
 - changed SCA_MS to always return an SCA with k1=k-1
 - made CMMSSY work as a stand-alone starter construction (always returns SCA with k1=k-1)
 - added further cover starters from LCDST paper, most to LCDST and the shorter ones to CMMSSY 
 - copied suitable f=1 length k cover starters of LCDST to CMMSSY, where they work for k+1
 - improved documentation for CS construction data
 - bug fix: productPCA did not work for PCAs for which the top right part had level>=i in row i
 - bug fix: CA_to_PCA removed the attributes from D for the cheap variant without tryhard=TRUE

August 30 2025, version 0.17
 - updated the data frame powerCTcat with 179 additional constructions, based 
   on new implementations in the recent past; also, one bug with a claim of
   a too small design was fixed; several runs are needed, as now some ingredients 
   are of type powerCT
 - fixed the creation file for powerCTcat in extdata
 - added nconst columns to ColbournKeriCombis and PALEYcat, and removed redundant rows from 
   ColbournKeriCombis
 - bug fix: CS_CK did not yield class ca with attributes

August 27 2025, version 0.16
 - renamed smcCA to scphfCA, and made its workhorse function internal; adapted
   everything to that change
 - modified the workhorse function SMC to avoid the creation of unnecessarily
   large matrices, and to handle both the SCPHFs from the Sherwood et al. paper
   (type="2006") and from Lanus and co-authors (type="2018") by taking care 
   of their respective GFs for v=8 and v=9
 - added the Lanus SCPHFs for up to v=9, as far as they improve constructions,
   in the list CL_SCPHFs
 - provided further internal functions for Galois field calculations, e.g., 
   the determinant, for checking CPHFs
 - fixed sloppy writing in the documentation of MTVTRouxtypeCA
   
August 21 2025, version 0.15
 - implemented one more column for strength 3 with v=2^s in SCA_Busht, 
   with consequences also for fuse_BushtCA, and fixed the manual for 
   SCA_Busht and SCA_Bose
 - adapted the N and k functions for SCA_Busht and fuse_BushtCA
 - implemented the Sherwood, Martirosyan, Colbourn (2006) construction
   with CPHFs and permutation vectors
 - created new file auxiliary_GFcalculus.R that holds all Galois field 
   related auxiliary functions (some moved from auxiliary_functions.R),
   including a somewhat hacky mygf for handling cases that require a 
   different characteristic polynomial than that used in lhs
 - improved the documentation of the verbose=2 case for function coverage
 
August 15 2025, version 0.14

 - added constructions by Cohen, Colbourn and Ling (2003, 2008) based 
   on ordered designs, in function ODbasedCA, with N_ and k_ functions,
   and made it available in Ns and ks and thus bestN and bestCA
 - added a CA(12,3,8,2) with two constant runs (available as ca12.2.8, 
   or via miscCA(3,8,2))
 - adapted the power construction info to the availability of the CA(12,3,8,2)
   with two constant rows
 - removed three arrays from TJ2level_CAs, and adapted corresponding rows of TJcat;
   added replacement constructions also for arrays that were never provided because 
   of their size, and removed N results for those
 - adapted N_TJcat and k_TJcat to yield result that refers to stored arrays and
   ignore removed better CAs from other constructions; this way, the actual 
   usage of stored arrays is more transparent
 - improved and renamed now internal Martirosyan and van Trung (2004) 
   Roux type implementations, and their referencing of actual theorems in the paper; 
   added a an API function MTVTRouxtypeCA that takes t, k, v as inputs, 
   and added the function N_upper_MTVTRouxtypeCA
   for the size (which can only yield upper bounds because of 
   unforeseeable number of duplicates that can be removed);
   not included in Ns and thus in bestCA so far (and probably rarely best, 
   if ever)
 - moved SCA_Busht further up in Ns, because its arrays yield better results
   than OAs in miscCA when used in recursive constructions (via bestCA)
 - modified miscCA to not make rows constant per default, 
   because this can make large arrays very slow
 - renamed internal function OD to OD2 and moved it to auxiliary_ingredients
 - bug fix: made bestCA obey to its argument fixNA
   
August 11, 2025, version 0.13

 - adapted projectionBose to also include the case with s=q+1 and c=1
 - created a function projBoseMCA for mixed level CAs
 - created a function N_upper_MCA
 - removed the futile run size warning for CS_LCDST with f>3
 - bug fix: N_CS_LCDST allows k smaller than min(k_implemented), 
   where k_implemented is the smallest k 
   for which there is a starter for the given v;
   CS_LCDST did not work for k < min(k_implemented); 
   fixed by making it work the same way as CS_MS
 - bug fix: k_LCDST did not work correctly (too small) for f > 3 
 - bug fix: CS_LCDST was not accessible by bestCA because of a typo


08 August 2025, version 0.12

 - added CS_CMMSSY for direct cover starters; the arrays have the same sizes as 
   CS_MS, but are guaranteed to be SCAs and can therefore be used recursively
   with function productPCA (similar to recursiveBose, also implemented in 
   function CS_CMMSSY)
 - added function fuse_BushtCA and corresponding N and k functions
 - added function compositCA and corresponding N and k functions 
   (for automated composition, using work horse function crossCAs)
 - modified N_powerCA and k_powerCA to allow t > requested t (will often yield 
   unattractive numbers; may be reverted after gaining experience)
 - added attributes to output of SCA_Busht
 - bug fix for crossCAs: construction did not have the correct number of repeats,
   and yielded a wrong result with a warning
 - bug fix for powerCA: cases with implemented t larger than requested t yielded an 
   uncaptured error
 - bug fix for dwyerCA: did not work for k less than maximum possible number of columns
   
04 Aug 2025, version 0.11
 - implemented an iterative version of function coverage
   that also allows to abort at the first non-perfect coverage
 - fixed a bug that was introduced on 31 July that made coverage fail
   for mixed level CAs

31 July 2025, still version 0.10
 - documented clearly that nlevels order does not determine the output 
   column order in MCA2
 - hand over objects from the function environment to the parallel cluster
   in function coverage
 - Bug fix: Ns_CK_doubling was wrong for v=3, as it did not multiply plus
   with v-1

23 July 2025, version 0.10
 - modified postopNCK to include an outer loop 
   (new arguments outerRetrymax and outermaxnochange, and seed),
   and deferred the previous inner workings to an internal function
 - provide detailed iteration information of postopNCK by messages (suppressible),
   stating the interim run size achieved at the end of each outer loop, 
   supporting informed decision on early interruption of the process
 - permit interruption of postopNCK with interim result preserved
 - implemented stop in postopNCK when lower bound on run size is reached
 - added function MCA2 for creating strength 2 mixed CAs with a generalization 
   of the first method of Sherwood (2008), as stated in Groemping (2025)
 - made CA_to_MCA interruptable with latest result returned
   (as the run size reduction may be slow)
 - implemented the possibility for outerRetry=0 (suppressing optimization)
 - made CA_to_MCA and MCA2 return objects with meaningful attributes; 
   handling of the information from the ingoing array may benefit from 
   future improvements
 - bugfix: postopNCK did not fix the top fixedrows rows; now it does
 - bugfix: CA_to_MCA was not reproducible for a given seed; now it is

13 July 2025, version 0.9
 - added function CA_to_MCA for creating mixed level CAs
 - changed the output of function postopNCK in case the
   number of rows cannot be reduced
 - fixed bug in function postopNCK that led to violation of 
   the requested strength in some cases

10 July 2025, version 0.8
 - implemented data file powerCTcat in support of power
   construction according to Colbourn and Torres-Jimenez (2010)
 - implemented functions powerCA (with type="CT"), N_powerCT and 
   k_powerCT and incorporated them into Ns, bestN, bestCA
 - implemented full factorial for k<=t (which sometimes occurs  
   for ingredients in recursive constructions)
 - omitted exclude="CK_doublingCA" for CK_doubling ingredients,
   which substantially increases the covered strength 3 scenarii
   (though not always with competitive arrays)
 - removed seven arrays from TJcat and TJ2level_CAs, as well 
   as from the folder TJ2level in extdata
 - moved SCA_Bose and SCA_Busht to a common file with dedicated
   documentation
 - fixed various bugs in other constructions that showed up in 
   testing the new function powerCA, among them the number 
   of columns for fuseBoseCA, which was always the maximum 
   instead of the requested k

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
