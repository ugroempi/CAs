### check the dimensions of the CAs created by scphfCA
for (i in 1:nrow(SCPHFcat)){
  print(i)
  D <- scphfCA(SCPHFcat$t[i], SCPHFcat$k[i], SCPHFcat$v[i])
  if (!nrow(D)==SCPHFcat$N[i] && ncol(D)==SCPHFcat$k[i]) stop("wrong dimensions")
  if (!max(D)==SCPHFcat$v[i]-1) stop("wrong v")
}
