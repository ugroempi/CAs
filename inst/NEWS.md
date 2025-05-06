# NEWS

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
