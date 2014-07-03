landgenreport <- function(cats,
                          fric.raster,  #friction matrix   
                          gen.distance = "Gst.Nei",  #Gst Hedrick, Gst_Nei, smouse, kosman, D, %allels shared 
                          NN=4, 
                          pathtype="leastcost",
                          plotpath=TRUE, 
                          theta=1, #for randomSP distance
                          mk.resistance=TRUE,    #mantel tests 
                        
                            mapdotcolor ="blue",
                            mapdotsize=1,
                            mapdotalpha=0.4,
                            mapdottype=19 ,
                            mapzoom=NULL,                          
# "roadmap","mobile","satellite","terrain","hybrid","mapmaker-roadmap","mapmaker-hybrid"

                          
                          mk.custom = FALSE,
                          

                          fname="LandGenReport",
                          foldername="results",
                          path.pgr=NULL,
                          mk.Rcode=FALSE,       # make the code that was ran available as an R file
                          mk.complete=FALSE,    # create a full report)  
                          mk.pdf=TRUE)
{
  if (class(cats)!="genind") {cat("You did not provide a valid genind object! Script stopped!\n"); return;}
  
  # Check for combinations of populations and loci with low numbers of individuals and alleles  
  npops<-length(levels(cats@pop))
  nloci<-length(cats@loc.names)
  
  # this splits bilby up into loci
  loci<-seploc(cats)
  
  # this further subdivides the loci into populations
  locipop<-lapply(loci,seppop)
  
  popsizes<-matrix(NA,nrow=nloci,ncol=npops)
  for (i in 1:nloci){
    for (j in 1:npops){
      popsizes[i,j]<-sum(!is.na(apply(locipop[[i]][[j]]@tab,1,sum)))
    }
  }
  
  for(i in 1:dim(popsizes)[2]){
    numlow<-length(which(popsizes[,i]<3))
    if(numlow>0) message("Population ",unname(cats@pop.names)[i]," has ",numlow," locus/loci with less than 3 genotypes. This may cause errors in some analyses. We advice to combine or drop populations with low numbers of genotypes. ")
  }
  
  
  #cut down length of loci names to  6  and make sure they are unique
  cats@loc.names <- substr(cats@loc.names,1,6)   
  if (length(unique(cats@loc.names))!= length(cats@loc.names)) 
  {cats@loc.names <- paste(1:length(cats@loc.names),"-",substr(cats@loc.names,1,4), sep="")
  cat("Loci names were not unique and therefore adjusted.\n")
  }

#check if ind.names are unique!!!!
#adjust if necessary and issue a notification
if (length(unique(cats@ind.names))!=length(cats@ind.names)) 
  {cats@ind.names <- paste(1:length(cats@ind.names),"-",substr(cats@ind.names,1,8),sep="")
  cat("Individual names were not unique and therefore adjusted.\n")
  }


#check if pop.names are unique!!!!
#adjust if necessary and issue a notification  
if (length(unique(cats@pop.names))!=length(cats@pop.names)) 
  {
  cats@pop.names <- paste(1:length(cats@pop.names),"-",substr(cats@pop.names,1,6),sep="")
  cat("Subpopulation names were not unique and therefore adjusted.\n")
  }

 
#set directory where to save a file, defaults to tempdir (follow R policy)
  if (is.null(path.pgr)) 
  {
  path.pgr <- tempdir()
  }

  
  #create a foldername folder if not existing...
  dirfiles <- list.dirs(path=path.pgr, recursive=FALSE)
  if (!(tolower (file.path(path.pgr,foldername))) %in% tolower(dirfiles)) {
    dir.create(file.path(path.pgr,foldername))
    cat("There is no ",foldername, " folder. I am trying to create it; \notherwise please create the folder manually. \n")
  }
  owd <-getwd()
  setwd(file.path(path.pgr, foldername))
 
  # coordinates must be in xy !!!!!
  coords=FALSE
  if (is.null(cats@other$xy)) coords=FALSE 
  if ( nrow(cats@other$xy) == length(cats@ind.names) ) coords=TRUE
 
  # give cats a filename that can be seen in the snw chunks
  cats@other$filename<- fname
  cats@other$foldername<-foldername
  cats@other$path <- path.pgr
  #determine the type of map
  if (coords) 
  {

  cats@other$mapdotcolor =mapdotcolor
  cats@other$mapdotsize=mapdotsize
  cats@other$mapdotalpha=mapdotalpha
  cats@other$mapdottype=mapdottype
  cats@other$mapzoom=mapzoom
  
  }  
  
###################################
##### create a new environment to run knitr in it
 pgr <- new.env(parent=.GlobalEnv)
 assign("cats",cats,envir=pgr)


###################################



  # save the data in a tempfile
 # save(cats, file=paste(foldername,"\\","tempcats.rdata",sep=""))
  
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
if (is.null(path)) {cat("Could not find snw files in the PopGenReport library folder. Please check if the package is installed correctly [e.g.  installed.packages() ]. \n"); return;}
  #for testing:
  #path <- "d:\\bernd\\R\\popgenreport\\inst\\swchunks\\"
  #path<- "C:\\Aaron files\\popgenreport098\\PopGenReport_0.98\\PopGenReport\\swchunks\\"
  header.file <- readLines(paste(path,"header.snw",sep=""))
  required<- readLines(paste(path,"required.snw",sep=""))
  compl<-c(header.file,required) 
  
  cat("Compiling report...\n")
  if (coords==FALSE) cat(" - No valid coordinates were provided. \n   Be aware you need to provide a coordinate (or NA) for each individual\n   and the coordinate heading in slot @other has to be 'latlong' or 'xy'.\n   All analyses will be skipped!\n") 
 
if ((mk.resistance==TRUE | mk.complete==TRUE)  & (coords & !is.null(fric.raster))) 
  {
    cat("- Landscape genetic analysis using resistance matrices...\n")  
 #   fr.raster<<-fric.raster
    assign("gen.dist",gen.distance,envir=pgr)
    assign("fr.raster",fric.raster,envir=pgr)
    assign("NN",NN,envir=pgr)
    assign("pathtype",pathtype,envir=pgr)
    assign("theta",theta,envir=pgr)
    assign("plotpath",plotpath,envir=pgr)
    
 #   gen.dist <<- gen.distance
    pmantel<-  readLines(paste(path,"pmantel.snw",sep=""))
    compl<-c(compl,pmantel)
  } 
  

if (mk.custom==TRUE){
  cat("- Run customised snw file, custom.snw ...\n")
  custom<-readLines(paste(path,"custom.snw",sep=""))
  compl<-c(compl,custom)
}



footer.file<-readLines(paste(path,"footer.snw",sep=""))  
compl<-c(compl,footer.file)

#compl <- c(header.file, required, loaddata, mapping, popheterozygosity, footer.file)


rnwfile <- paste(fname,".rnw",sep="")
texfile <-  paste(fname,".tex",sep="") 


zz <- file(file.path(path.pgr,foldername,rnwfile), "w")
writeLines(compl,zz)
close(zz) 


#setwd(paste(path.pgr,foldername, sep="/"))
cat(paste("Analysing data ...\n", sep=""))
#Sweave(paste(fname,".rnw",sep=""), output=paste(fname,".tex",sep=""), quiet=FALSE, driver=mydriver)
flush.console()
knit(input=rnwfile, output=texfile, quiet=TRUE, envir=pgr)

if (mk.pdf==TRUE)
{
cat(paste("Creating pdf from: ",rnwfile," ...\n",sep=""))
knit2pdf(texfile, texfile)
cat(paste("Finished.\nCheck ",fname,".pdf for results.\n", sep=""))
}

if (mk.Rcode) {
  cat(paste("Creating R code from: ",rnwfile,"...\n"), sep="")
  rfile <-paste(fname,".R",sep="")
  purl(input=rnwfile, output=rfile)
#  Stangle(paste(fname,".rnw",sep=""))
}
    

cat(paste("All files are available in the folder: \n",file.path(path.pgr, foldername),"\n",sep=""))

#reset working directory to previous
setwd(owd)
return(pgr$allresults)
}
