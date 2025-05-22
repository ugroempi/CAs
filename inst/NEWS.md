# NEWS

22 May 2025, version 0.2.2

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
