gd_kosman<-function (population, verbose = TRUE) 
{
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }
  maxind <- (dim(population@tab))[1]
  maxloci <- length(population@loc.names)
  d.fast <- matrix(NA, nrow = maxind, ncol = maxind)
  loci.used <- matrix(NA, nrow = maxind, ncol = maxind)
  cs <- cumsum(population@loc.nall)
  cs.l <- c(1, cs[-length(cs)] + 1)
  for (i in 1:(maxind - 1)) {
    if (verbose) 
      message("i=", i)
    for (j in i:maxind) {
      if (i == j) {
        d.fast[j, i] <- 0
      }
      else {
        i1 <- population@tab[i, ]
        i2 <- population@tab[j, ]
        nas <- sum(is.na(sapply(1:maxloci, function(x, 
                                                    cs, cs.l, i1, i2) sum(i1[cs.l[x]:cs[x]] + i2[cs.l[x]:cs[x]]), 
                                cs, cs.l, i1, i2)))
        d.fast[j, i] <- 1 - (sum(i1 == 0.5 & i2 == 0.5, 
                                 na.rm = T) * 0.5 + sum(i1 == 1 & i2 == 1, na.rm = T) + 
                               sum(i1 == 1 & i2 == 0.5, na.rm = T) * 0.5 + 
                               sum(i1 == 0.5 & i2 == 1, na.rm = T) * 0.5)/(maxloci - nas)
        loci.used[j,i]<-maxloci-nas
      }
    }
  }
  colnames(d.fast) <- population@ind.names
  rownames(d.fast) <- population@ind.names
  d.fast[upper.tri(d.fast, diag = FALSE)] <- NA
  #d.fast <- as.dist(d.fast)
  #loci.used<-as.dist(loci.used)
  kosman.out<-list(geneticdist=d.fast,loci_used=loci.used)
  return(kosman.out)
}