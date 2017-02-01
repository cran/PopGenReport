#' Function to convert textfiles into a \linkS4class{genind} object (the format
#' required for popgenreport)
#' 
#' This function converts genetic data provided in a tabular 'comma separated
#' value' (.csv) file into a genind object, the data format required to run
#' PopGenReport. At the moment it works only for codominant markers (e.g.
#' microsatellites). This function is based on df2genind from the
#' \code{\link{adegenet}} package.
#' 
#' The format of the .csv file is very important. Make sure your \bold{headings
#' are exactly as provided in the example file} or the conversion will likely
#' fail (e.g. use lower cases for all headings). Use your favourite text editor
#' to reformat the file (e.g. Excel) to prepare the data and save it as csv
#' file. You need to provide the number of columns for each of your data
#' sections. These are: ind, pop, lat, long, other.min, other.max, and whether
#' there is a single column per allele (if you use a single column for both
#' alleles then you need to specify the seperator as well), or two columns per
#' allele. Please refer to the example files to make sure your file is in the
#' correct format and then check the conversion by typing:
#' 
#' mydata <- read.genetable(\"mygeneticdat.csv\") mydata
#' 
#' The easiest way to provide spatial coordinates is to use the read.genetable
#' function and use the \code{lat}, \code{long} arguments for WGS1984 projected
#' data (the most common projection globally). For additional information how
#' to use spatial data with PopGenReport refer to the help of the
#' \code{popgenreport} function and to the popgenreport manual. \cr === There
#' are 3 treatments for missing values === - NA: kept as NA.
#' 
#' - 0: allelic frequencies are set to 0 on all alleles of the concerned locus.
#' Recommended for a PCA on compositionnal data.
#' 
#' - \"mean\": missing values are replaced by the mean frequency of the
#' corresponding allele, computed on the whole set of individuals. Recommended
#' for a centred PCA.
#' 
#' === Details for the sep argument === this character is directly used in
#' reguar expressions like gsub, and thus require some characters to be
#' preceeded by double backslashes. For instance, \"/\" works but \"|\" must be
#' coded as \"\\|\".
#' 
#' @param filename the name of your file in .csv file format
#' @param pop the column number if subpopulations are known, otherwise set to NULL
#' @param ind the column number of an individual identifier, otherwise set to NULL
#' @param lat the column number where the latitudinal coordinate is recorded
#' (can be set to NULL)
#' @param long the column number where the longitudinal coordinate is recorded
#' (can be set to NULL)
#' @param x the column number where the x coordinate is recorded. If not in
#' Mercator it needs to be transformed!!! (can be set to NULL)
#' @param y the column number where the y coordinate is recorded. If not in
#' Mercator it needs to be transformed!!! (can be set to NULL)
#' @param other.min if your data has some additional data (e.g. gender, size
#' etc.) you can include this data in the genind object as well. Values in this
#' section have to be in adjacent columns with no other columns between them.
#' other.min is the column number where this section starts.
#' @param other.max if your data has some additional data (e.g. gender, size
#' etc.) then you can convert this data as well. This section has to be in a
#' consecutive order with no other type of columns in between. other.max is the
#' column number where this section ends.
#' @param oneColPerAll needs to be specified. If your data is coded as a single
#' column per allele
#' @param NA.char can be NA, 0 or mean. See details section.
#' @param sep a character string separating alleles. See details.
#' @param ncode an optional integer giving the number of characters used for coding one genotype at one locus.
#'  If not provided, this is determined from data.
#' @param ploidy Ploidy of the data set. Be aware that most analysis do only
#' run for diploid data and that missing data in polyploid data sets are
#' ambigious, which may give dubious results if not handled appropriately.
#' @return an object of the class genind This kind of object is needed to be
#' passed on to the popgen.report function.
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au)
#' @seealso \code{\link{import2genind}}, \code{\link{df2genind}},
#' \code{\link{read.fstat}}, \code{\link{read.structure}},
#' \code{\link{read.genetix}} \code{\link{read.genepop}}
#' @examples
#' 
#' #example file with one column per loci, seperated by forwardslash
#' read.csv(paste(.libPaths()[1],"/PopGenReport/extdata/platypus1c.csv", sep="" ))
#' platy1c <- read.genetable( paste(.libPaths()[1],"/PopGenReport/extdata/platypus1c.csv"
#' , sep="" ), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=FALSE, 
#' sep="/", )
#' 
#' 
#' #example file with two columns per loci
#' read.csv(paste(.libPaths()[1],"/PopGenReport/extdata/platypus2c.csv", sep="" ))
#' platy2c <- read.genetable( paste(.libPaths()[1],"/PopGenReport/extdata/platypus2c.csv",
#'  sep="" ), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=TRUE)
#' 
#' #to get a pdf output you need to have a running Latex version installed on your system.
#' #run a report (with a map)
#' #res<- popgenreport(platy2c, mk.counts=TRUE, mk.map=TRUE, mk.allele.dist=TRUE, mk.pdf=TRUE)
#' @export
#' 
read.genetable <- function(filename, pop=NULL, ind=NULL,lat=NULL, long=NULL, x=NULL, y=NULL, other.min=NULL, other.max=NULL,oneColPerAll,NA.char=NA,  sep=NULL, ncode=NULL ,ploidy=2 )
{
gfile <-read.csv(filename)

popnr=NULL
pops=NULL
indnr=NULL
inds=NULL

latlong=NULL

if ("pop" %in% colnames(gfile))            
  {
  pops <- as.factor(gfile$pop)
  popnr <- which(colnames(gfile)==tolower("pop"))
  }
if ("ind" %in% colnames(gfile)) 
  {
  inds <- as.factor(gfile$ind)
  indnr <- which(colnames(gfile)==tolower("ind"))
  }





genes <- gfile
rem=NULL

if (!is.null(other.min) & !is.null(other.max)) rem <- other.min:other.max         
if (!is.null(popnr)) rem <- c(rem, popnr)
if (!is.null(indnr)) rem <- c(rem, indnr)
if (!is.null(lat))  rem <- c(rem, lat)
if (!is.null(long)) rem <- c(rem, long)
if (!is.null(x)) rem <- c(rem, x)
if (!is.null(y)) rem <- c(rem, y)


if (!is.null(rem)) genes <- gfile[, -c(rem)]   else genes <- gfile


if (oneColPerAll==TRUE)
{
res <-data.frame(allele= rep(NA,dim(genes)[1]) )

for (i in seq(1,dim(genes)[2]-1,ploidy))
  {
  dummy <-genes[,i]
  for (ii in (i+1):(i+ploidy-1))
  {
  dummy<- paste(dummy,genes[,ii], sep="/")
  
  }
  
  res[,ceiling(i/ploidy)] <-  dummy
  
 

    colnames(res)[ceiling(i/ploidy)] <-  colnames(genes)[i]
  }

  
sep="/"

genes <- res
} 



  
  
  df <- df2genind(genes, pop=as.character(pops), ind.names=as.character(inds), NA.char=NA.char,     ncode=ncode,loc.names=colnames(genes), sep=sep, ploidy=ploidy)
  
  if (!is.null(lat) & !is.null(long))
  df@other$latlong <- data.frame(lat=gfile[,lat], long=gfile[,long])
  if (!is.null(x) & !is.null(y))
  df@other$xy <- data.frame(x=gfile[,x], y=gfile[,y])
  
  
  if (!is.null(other.min) & !is.null(other.max)) 
  {
  df@other$data  <- data.frame(gfile[,c(other.min:other.max)])
  }

return(df)



}
