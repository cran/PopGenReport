gd_kosman <- function(population, verbose=TRUE){
  
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }
  # what is the total number of individuals
  maxind<-(dim(population@tab))[1]
  
  # what is the maximum number of loci
  maxloci<-length(population@loc.names)
  
  # build a full, symmetric dissimilarity matrix to fill
  d.fast<-matrix(NA,nrow=maxind,ncol=maxind)
  cs <- cumsum(population@loc.nall) # this determines the end of each loci frame for each ind's genotype
  cs.l <- c(1,cs[-length(cs)]+1) # this determine the beginning of each loci frame for each ind's genotype
  
  for (i in 1:(maxind-1)){
    if (verbose) message("i=",i)
    for (j in i:maxind){ 
      if (i==j){      # an individual isn't dissimilar from itself
        d.fast[j,i]=0
      }
      else {      # comparing two different individuals 
        i1 <- population@tab[i,]
        i2 <- population@tab[j,]
        
        # Kosman and Leonard 2005
        nas <- sum(is.na(sapply(1:maxloci, function(x,cs,cs.l,i1,i2) sum(i1[cs.l[x]:cs[x]] +i2[cs.l[x]:cs[x]])  , cs, cs.l,i1,i2)))
        d.fast[j,i] <- 1-(sum(i1==0.5 & i2==0.5, na.rm=T)*0.5 + sum(i1==1 & i2==1, na.rm=T) +   sum(i1==1 & i2==0.5, na.rm=T)*0.5 + sum(i1==0.5 &i2==1, na.rm=T)*0.5) / (maxloci-nas)
      }
    }
  }
  
  # put names on the rows and columns on d.fast
  colnames(d.fast)<-population@ind.names
  rownames(d.fast)<-population@ind.names
  
  # force upper triangle to NA
  d.fast[upper.tri(d.fast,diag=FALSE)]<-NA
  d.fast<-as.dist(d.fast)
  return(d.fast)
}