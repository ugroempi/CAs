## create_rmsCOPIEScat
  ## Combinations (q, m) whose arrays, as produced by RMScopies(q, m), were
  ## verified exhaustively (October 2026):
  ## - all except q = 19: the R output was checked with ca_verify.c
  ##   (C kernel of caverify, from github.com/chrispbsmolen/ca-tools);
  ## - q = 19: the R output coincides (same multiset of rows) with the arrays
  ##   of the Python implementation, which passed an independent brute-force
  ##   checker (verify.c in github.com/MichelBacon/covering-arrays-rms-extension).
  aus <- rbind(
    data.frame(q = 3,  m = 0:13),
    data.frame(q = 4,  m = 0:21),
    data.frame(q = 5,  m = 0:31),
    data.frame(q = 7,  m = 0:57),
    data.frame(q = 8,  m = 0:73),
    data.frame(q = 9,  m = 0:91),
    data.frame(q = 11, m = c(0:13, 22)),
    data.frame(q = 13, m = 0:14),
    data.frame(q = 16, m = c(0, 1, 16)),
    data.frame(q = 19, m = c(1, 19, 20))
  )
  aus$k <- aus$q^2 + aus$q + 1 + aus$m
  aus$N <- mapply(.RMS_N, aus$q, aus$m)
  aus$verifier <- ifelse(aus$q == 19, "identical to verified Python array",
                         "ca_verify (ca-tools)")
rmsCOPIEScat <- aus
save(rmsCOPIEScat, file="data/rmsCOPIEScat.rda", compress="xz")

