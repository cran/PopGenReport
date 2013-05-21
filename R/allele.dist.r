allele.dist<-function(population, mk.figures=TRUE){
  # package require adegenet and pegas
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }

  # initial steps...
  numloci<-length(population@loc.nall)  # this gets the total number of loci across all pops
  numpops<-length(population@pop.names) # this gets the total number of pops
  popnumallele<-population@loc.nall     # this is a list of the population wide number of alleles at each pop
  lociname<-attributes(popnumallele)[[1]] # this is a list of the locinames (just L01, L02, L03,...)
  subdivpops<-seppop(population)

  # create list of matrices in which to place the numbers from summary
  alleletable<-vector("list",numloci)
  fralleletable<-vector("list",numloci)
  for(i in 1:numloci){
    alleletable[[i]]<-matrix(nrow=popnumallele[[i]],ncol=numpops)
    colnames(alleletable[[i]])<-population@pop.names
    rownames(alleletable[[i]])<-population@all.names[[i]]
    fralleletable[[i]]<-matrix(nrow=popnumallele[[i]],ncol=numpops)
    colnames(fralleletable[[i]])<-population@pop.names
    rownames(fralleletable[[i]])<-population@all.names[[i]]
  }

  # this is going to loop over all populations
  for (i in 1:numpops){
    x<-as.loci(subdivpops[[i]])
    s<-summary(x)
    # this loops over the loci
    for (j in 1:numloci){
    # this is the number of 
      namevec<-as.numeric(names(s[[j]]$allele))
      numnames<-length(namevec)
      # j<-2
      tablenames<-as.numeric(rownames(alleletable[[j]]))
      for (k in 1:numnames){
        rownum<-which(tablenames==namevec[k])
        #  message("i = ",i," j = ",j," k = ",k," rownum = ",rownum)
        alleletable[[j]][rownum,i]<-s[[j]]$allele[k]
      }  
    }
  }

  allpops<-as.loci(population)
  numbers<-summary(allpops)
  checkcnts<-matrix(nrow=numloci,ncol=2)
  for (i in 1:numloci){
    checkcnts[i,1]<-sum(numbers[[i]]$allele)
    checkcnts[i,2]<-sum(alleletable[[i]],na.rm=T)
  }

  for (i in 1:numloci){
    for (j in 1:numpops){
      colsum<-sum(alleletable[[i]][,j],na.rm=T)
      fralleletable[[i]][,j]<-round(alleletable[[i]][,j]/colsum, digits=3)
    }
  }
  
  if (mk.figures){
    breaks<-seq(0,1,.05)
    color.palette  <- colorRampPalette(c("yellow", "red"))(length(breaks) - 1)
    for (i in 1:numloci){
      figlabel<-paste("Loci: ",population@loc.names[i]," List # ",i,sep="")
      heatmap.2(fralleletable[[i]], dendrogram="none", Rowv=NA, Colv=NA,col = color.palette,breaks=breaks, scale="none", trace="none", cellnote=alleletable[[i]], margins=c(5,10),main=figlabel, lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 4, 2 ), notecol="black",tracecol="black",linecol="black")
    }
  }
  alleletables<-list(count=alleletable, frequency=fralleletable)
  return(alleletables)
}