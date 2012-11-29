popgenreport <- function(cats=NULL,
                          
                          mk.counts=T,   # this switch is to provide a population overview
                          mk.map=F,        # this switch is for the map
                            maptype="satellite",
                            mapdotcolor ="blue",
                            mapdotsize=1,
                            mapdotalpha=0.4,
                            mapdottype=19 ,
                          
#c("roadmap",
#"mobile",
#"satellite",
#"terrain",
#"hybrid",
#"mapmaker-roadmap",
#"mapmaker-hybrid")",
#                          
                          
                          mk.locihz=F,     # this switch is to test for population heterozygosity
                          mk.hwe=F,   # this switch is for population wide HWE
                          mk.subgroups=NULL,  # this switch is to run the analysis for subgroups
                          mk.fst=F,        # this switch is to run FST tests on the full population
                          mk.gd.smouse=F,   # this switch is to run the Smouse and Peakall genetic distances
                          mk.gd.kosman=F,   # this switch is to run the Kosman and Leonard genetic distances

                          mk.allele.dist=F, # this switch it to look at allele distributions by loci and pop

                          mk.differ.stats=F ,     # this switch is to look at population differentiation statistics (Fst, Gst, etc)
                          fname="PopGenReport",
                          path.pgr=NULL,
                          mk.Rcode=F,       # make the code that was ran available as an R file
                          mk.complete=F,    # create a full report)  
                          mk.pdf=T)
{
  if (class(cats)!="genind") {cat("You did not provide a valid genind object! Script stopped!\n"); return;}
  
    #set directory where to save a file, defaults to tempdir (follow R policy)
  current.dir <- getwd()
  if (is.null(path.pgr)) path.pgr <- tempdir()
  
  
  setwd(path.pgr)
  
  #create a figures folder if not existing...
  dirfiles <- list.dirs(recursive=F)
  if (!("./figures" %in% tolower(dirfiles))) {
    dir.create("figures")
    cat("There is no figures folder. I am trying to create it otherwise create it manually. \n")
  }

  # conversion of lat longs to google map data (Mercator (dismo) wants to have long lat)
  if (!is.null(cats@other$latlong)) cats@other$mercat <- Mercator(cats@other$latlong[,c(2,1)])
  
  # give cats a filename that can be seen in the snw chunks
  cats@other$filename<-fname
  
  #determine the type of map
  if (mk.map==T | mk.complete) 
  {
  cats@other$maptype=maptype
  cats@other$mapdotcolor =mapdotcolor
  cats@other$mapdotsize=mapdotsize
  cats@other$mapdotalpha=mapdotalpha
  cats@other$mapdottype=mapdottype
  }  
  
  #stick in group for subgroup routines
#if mk.subgroups==T use from group cats@other$data$group
if (!is.null(mk.subgroups))
{
if (mk.subgroups[1]==T & !is.null(cats$other$data$group)  )
 
{
cats@other$group <- cats$other$data$group
}
#if subgroups is a factor of the according length then use this...
 
if (is.factor(mk.subgroups))
{
cats@other$group <- mk.subgroups
mk.subgroups <-T
}
if (!is.null(mk.subgroups) & is.null(cats$other$group)) 
  {
  cat("A report for subgroups was requested (mk.subgroups was not NULL), but I could not find a subgroup definition. Set yourgenindobject@other$group to a appropriate factor. mk.subgroups set to FALSE)")
  mk.subgroups=F
  }
 
}  else mk.subgroups=F

  # save the data in a tempfile
  
  save(cats, file="tempcats.rdata")
  
  #check path to the snw files
path <- NULL
  for(i in seq_along(.libPaths()))
{
  if (file.exists(paste(.libPaths()[i],"/PopGenReport/swchunks/header.snw",sep="")))  
  {
  path <-   paste(.libPaths()[i],"/PopGenReport/swchunks/", sep="" )
  break
  }
  
}
if (is.null(path)) {cat("Could not find snw files in the PopGenReport library folder. Please check if the package is installed correctly (e.g.  installed.packages()[\"PopGenReport\",2]). \n"); return;}
  #for testing:
  #path <- "d:\\bernd\\R\\popgenreport\\inst\\swchunks\\"
  
  header.file <- readLines(paste(path,"header.snw",sep=""))
  required<- readLines(paste(path,"required.snw",sep=""))
  loaddata<- readLines(paste(path,"load.data.snw",sep=""))
  compl<-c(header.file,required,loaddata) # deleted dataprep from the last slot in this call
  
  
  cat("Compiling report...\n")
  if(mk.counts | mk.complete){
    cat("- Generall summary...\n")
    overview<-readLines(paste(path,"counts.snw",sep=""))
    compl<-c(compl,overview)
  }
  
  if (mk.map | mk.complete){
    cat("- Map of individuals...\n")  
    mapping<-  readLines(paste(path,"map.snw",sep=""))
    compl<-c(compl,mapping)
  }
  
  if (mk.locihz | mk.complete){
    cat("- Statistics on population heterogeneity ...\n")  
    popheterozygosity <- readLines(paste(path,"locihz.snw",sep=""))
    compl<-c(compl,popheterozygosity)
  }
  
  if (mk.allele.dist | mk.complete){
    cat("- Allelic distances ...\n")  
    numloci<-length(cats@loc.nall)
    alleledistn <- readLines(paste(path,"allele.dist.snw",sep=""))
    compl<-c(compl,alleledistn)
  for (i in 1:numloci){
    figfilelabel<-paste(path,"allele",i,".pdf",sep="")
    line1ad<-"\\begin{figure}[hptb]"
    line2ad<-"\n"
    line3ad<-paste("\\centerline {\\includegraphics[width=6in,height=4.5in] {figures/",fname,"-allele",i,".png}}",sep="")
    line4ad<-"\n"
    line5ad<-"\\end{figure}"
    line6ad<-"\n"
    if (i%%10==0){
      line7ad<-"\\FloatBarrier"
      line8ad<-"\n"
      compl<-c(compl,line1ad,line2ad,line3ad,line4ad,line5ad,line6ad,line7ad,line8ad)
    } else {
      compl<-c(compl,line1ad,line2ad,line3ad,line4ad,line5ad,line6ad)
    }
    compl<-c(compl,"\n")
  }  
}

if (mk.fst| mk.complete){
     cat("- Pairwise Fst ...\n")  
  popfst<-readLines(paste(path,"fst.snw",sep=""))
  compl<-c(compl,popfst)
}
if (mk.differ.stats | mk.complete){
     cat("- Pairwise differentiations ...\n")  
  differentiate<-readLines(paste(path,"differ.stats.snw",sep=""))
  compl<-c(compl,differentiate)
}

if (mk.hwe | mk.complete){
  cat("- Test for Hardy-Weinberg-Equilibria ...\n")  
  popHWEll<-readLines(paste(path,"hwe.snw",sep=""))
  compl<-c(compl,popHWEll)
}

if (mk.gd.kosman | mk.complete){
  cat("- Kosman & Leonard genetic distances...\n")
  kosman<-readLines(paste(path,"gd.kosman.snw",sep=""))
  compl<-c(compl,kosman)
}

if (mk.gd.smouse | mk.complete){
  cat("- Smouse & Peakall genetic distances...\n")
  smouse<-readLines(paste(path,"gd.smouse.snw",sep=""))
  compl<-c(compl,smouse)
}

if (mk.subgroups){
  cat("- Analysis for user defined subgroups...\n") 
  subgrp_header <- readLines(paste(path,"subgroups.header.snw",sep=""))
  compl<-c(compl,subgrp_header)
  # first get a list of the number of sexes
  sexlist<-attributes(cats@other$group)$levels
  # now count the number of distinct sexes tracked
  numsexes<-length(sexlist)
  subsec<-rep(NA,numsexes)
  for (j in 1:numsexes){
    i<-j           
    line1rpt<-"<<echo=false,results=hide>>="
    line2rpt<-"\n"
    line3rpt<-paste("i<-",i,sep="")
    line4rpt<-"\n"
    line5rpt<-"@"
    line6rpt<-"\n"
    compl<-c(compl,line1rpt,line2rpt,line3rpt,line4rpt,line5rpt,line6rpt)
    subsectcode<-readLines(paste(path,"subgroups.sectionout.snw",sep=""))
    compl<-c(compl,subsectcode)
  }
}
 

cat("\n") 
footer.file<-readLines(paste(path,"footer.snw",sep=""))
compl<-c(compl,footer.file)

#compl <- c(header.file, required, loaddata, mapping, popheterozygosity, footer.file)




zz <- file(paste(fname,".rnw",sep=""), "w")
writeLines(compl,zz)
close(zz) 

cat(paste("Analysing data ...\n", sep=""))
Sweave(paste(fname,".rnw",sep=""), output=paste(fname,".tex",sep=""), quiet=T)

if (mk.pdf==T)
{
cat(paste("Creating pdf from: ",fname,".rnw ...\n",sep=""))
tools::texi2pdf(paste(fname,".tex",sep=""))
cat(paste("Finished.\nCheck ",fname,".pdf for results.\n", sep=""))
}
if (mk.Rcode) {
  cat(paste("Creating Rcode: ", fname,".rnw ...\n"), sep="")
  Stangle(paste(fname,".rnw",sep=""))
}
    


cat("All figures are available in the figures subfolder\n")

cat(paste("Files are saved under: ", path.pgr,"\n", sep=""))

#reset working directory to previous

setwd(current.dir)

}
 
 
                                              