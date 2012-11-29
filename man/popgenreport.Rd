\name{popgenreport}
\alias{popgenreport}
\title{
This is the main function of the package. It analyses an object of class \code{\link{genind}} and then creates a report containing the results of the analysis. There are several routines that can be optionally included in the analysis and there are multip output options including a PDF with the report, R-code and an object (allresults) containing all of the results, which can be used for further analyses. 
}
\description{
This function is used to analyse population genetic data. The main idea is to provide a framework for analysing microsatellite genetic data using a mix of existing and new functions. The function works on an object of class \code{genind}. There are several ways to convert data into a \code{\link{genind}} object using existing functions provided by the \code{adegenet} package (  \code{\link{import2genind}}, \code{\link{df2genind}},\code{\link{read.fstat}}, \code{\link{read.structure}},
  \code{\link{read.genetix}} ,\code{\link{read.genepop}}) or refer to \code{\link{read.genetable}}  how to import data from an EXCEL (csv) document. The function performs a number of different genetic analyses (e.g. counts of indivuals and alleles across sub-populations, tests for heterozygosity and Hardy-Weinberg Equilibrium, differentiation statistics Fst, G'st, Jost's D, and genetic distance between individuals and populations), with users having the option to select which analysis routines are included in the report. To select a routine, the user simply turns on a switch e.g. mk.map=TRUE returns a map with the sampling location for each individual (if coordinates (latitudes and longitudes in WGS1984) are provided). 
}
\usage{
popgenreport(cats = NULL, mk.counts = T, mk.map=F, maptype="satellite", 
mapdotcolor ="blue", mapdotsize=1, mapdotalpha=0.4, mapdottype=19 , 
mk.locihz = F, mk.hwe = F, mk.subgroups = NULL, mk.fst = F, mk.gd.smouse = F, 
mk.gd.kosman = F, mk.allele.dist = F, mk.differ.stats=F, 
fname = "PopGenReport", path.pgr=NULL, mk.Rcode = F, mk.complete=F, mk.pdf=T)
}
\arguments{
  \item{cats}{this is the \code{\link{genind}} object the analysis will be based on. }
  \item{mk.counts}{
switch is to provide overview counts of the  number of individuals sampled, numbers of individuals and alleles sampled per sub-population, number of alleles per locus, mean number of alleles per locus and the percentatge of missing data.
}
  \item{mk.map}{switch to produce a map with the sampling location of each individual marked. This switch requires individual coordinates (latitudes and longitudes in WGS1984) be provided (under cats@other$latlong or see \code{\link{read.genetable}} on how to import them from a table of genetic data). An error message will be generated if you turn this routine on, but do not provide the coordinates in the right format. If the coordinates are provided in a seperate file, they must be attached to the genind object in the slot \cr\code{yourgenindobject@other$latlong <- yourlatlongdata}.\cr \code{yourlatlongdata} needs to be a data frame that has the same number and order of individuals per row as the population genetic data. Note that an internet connection is required to connect to the Google Maps server which provides the basemap for this routine.}
  \item{maptype}{Defines the type of map. Default is 'satellite'. Other options are: 'roadmap', 'mobile', 'terrain', 'hybrid'. }
  \item{mapdotcolor}{Color of dots for each individual on the map. Default is 'blue'.}
  \item{mapdotsize}{Size of dots for each individual. Default is 1.}
  \item{mapdotalpha}{Transparency of dots. 1 is invisible, 0 is no transparency. Default is 0.4.}
  \item{mapdottype}{Defines the type of the symbol. For explanation see pch under \code{\link{par}}. Default is 19 - a filled circle.}
  \item{mk.locihz}{switch to test for population heterozygosity}
  \item{mk.hwe}{switch to test for Hardy-Weinberg equilibrium for each loci and population}
  \item{mk.subgroups}{switch to repeat the entire analysis for subgroups. Either provide a subgroup as a factor (needs to have the same number of rows as your genetic data set, or you can provide it via the slot \code{yourgenindobject@other$group} and set \code{rptsubgroup=T}}
  \item{mk.fst}{switch to calculate Fst values for each loci and pairwise Fst (Nei's 1973) over subpopulations}
  \item{mk.gd.smouse}{Individual pairwise genetic distances based on Smouse and Peakall (1999). Refer to \code{\link{gd_smouse}}}
  \item{mk.gd.kosman}{Individual pairwise genetic distances based on Kosman & Leonhard (2005). Refer to \code{\link{gd_kosman}}}
  \item{mk.allele.dist}{switch to look at allele distributions by loci and subpopulation}
  \item{mk.differ.stats}{switch to look at population differentiation statistics (Nei's Gst, Hedrick's Gst, and Jost's D)}
  \item{fname}{filename for the output files. Defauts to PopGenReport. Note that using a filename which includes a space in the name will result in the filename for each figure being printed out in the PDF report for each figure. Replacing the space with an underscore should prevent this from happening. }
  \item{path.pgr}{Folder where the output files are stored. Defaults to the temporary directory (\code{tempdir()}). If you want to store the output in another directory, simply provide the path here. e.g. \code{path.pgr=getwd()} saves it in your current working directory.}

  \item{mk.Rcode}{switch to get the full R script that is used to generate the report. A great way to get a very detailed insight on the kind of analysis and also an easy way to generate a script which you can customize for your analytical needs.}
  \item{mk.complete}{switch to create a full report using all of the routines (all switches are set to \code{TRUE}, except \code{mk.subgroups}).}
  \item{mk.pdf}{  switch to create a shiny pdf output. You need a working \bold{latex} version running on your system (e.g. MikTex (Windows) or Texmaker (Linux, MacOSX). For more information refer to the \code{\link{Sweave}} package and its manual.}
}
\value{
The function returns an object that has all of the results produced by this function in it. The structure of the object can be accessed via \code{str(allresults)}. The main slots in this object (if you ran a full report) are:\cr
\code{dataoverview,PopHet,Alleledist,Fst,HsHtdifferentiate,HWEresults,subgroups,GDKosman,GDSmouse}

Additional ouput is provided in the form of a PDF which will be saved to your current working directory, and maps and figures which will be placed in a figures subfolder which may be generated automatically in your current working directory. 
If \code{mk.Rcode=T} is set, an R file named yourgeninobject.R will be saved to your current working directory.  

}
\references{
Kosman E., Leonard K.J. 2005. Similarity coefficients for molecular markers in studies of genetic relationships between individuals for haploid, diploid, and polyploidy species. Molecular Ecology 14:415-424

Peakall R., Smouse P. 2012. GenAlEx 6.5: Genetic analysis in Excel. Population genetic software for teaching and research - an update. Bioinformatics 28:2537-2539

%%Adamack & Gruber (2012)
}
\author{
Aaron Adamack & Bernd Gruber, aaron.adamack@canberra.edu.au, bernd.gruber@canberra.edu.au
}
\seealso{
  \code{\link{adegenet}}, \code{\link{pegas}},
  \code{\link{mmod}}
  }
\examples{
data(bilby)
popgenreport(bilby[1:100], mk.counts=TRUE, mk.map=TRUE, mk.pdf=FALSE)

### RUN ONLY with a working Latex version installed, see mk.pdf
#popgenreport(bilby, mk.counts=TRUE, mk.map=TRUE, mk.pdf=TRUE)

# for a full report in a single pdf set mk.complete to TRUE
# popgenreport(bilby, mk.complete=TRUE)
}
