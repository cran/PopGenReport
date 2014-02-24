\name{landgenreport}
\alias{landgenreport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a landscape genetic report
}
\description{
This function is the landscape genetic version of the \code{\link{popgenreport}} function. It needs to be provided with a genind object with spatial coordinates, a friction map (raster) and a specification which type of genetic distance should be used.  
Once all three type of input are provided with the necessary input, a landscape genetic analysis using least cost path analysis is computed (see Cushman et al. 2010, Landguth et al. 2010). Depending on the genetic distance meassurement this is done on a subpopulation basis (D, Gst.Hedrick, Gst.Nei=Fst) or on an individual basis (Kosman, Smouse).


}
\usage{
landgenreport(cats = NULL, mapdotcolor = "blue", mapdotsize = 1, mapdotalpha = 0.4,
mapdottype = 19, mapzoom = NULL, mk.pmantel = FALSE, fric.raster = NULL, 
gen.distance = "Gst.Nei", mk.custom = FALSE, fname = "LandGenReport", foldername = 
"results", path.pgr = NULL, mk.Rcode = FALSE, mk.complete = FALSE, mk.pdf = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cats}{ a \code{genind} object with spatial coordinates in the other slot

}
  \item{mapdotcolor}{ see \code{\link{popgenreport}}
}
  \item{mapdotsize}{ see \code{\link{popgenreport}}
}
  \item{mapdotalpha}{ 
see\code{\link{popgenreport}}
}
  \item{mapdottype}{ see \code{\link{popgenreport}}
}
  \item{mapzoom}{ see \code{\link{popgenreport}}
}
  \item{mk.pmantel}{ switch to do the analysis, should be set to TRUE
}
  \item{fric.raster}{ friction raster, that specifies the landscape where the analysis should be computed on
}
  \item{gen.distance}{type of genetic distance that should be used. Depending on the genetic distance meassurement this is done on a subpopulation basis (D, Gst.Hedrick, Gst.Nei=Fst) or on an individual basis (Kosman, Smouse).
}
  \item{mk.custom}{ switch to add a customised part to the landgenreport
}
  \item{fname}{ see \code{\link{popgenreport}}
}
  \item{foldername}{ see \code{\link{popgenreport}}
}
  \item{path.pgr}{ see \code{\link{popgenreport}}
}
  \item{mk.Rcode}{ see \code{\link{popgenreport}}
}
  \item{mk.complete}{ see \code{\link{popgenreport}}
}
  \item{mk.pdf}{see \code{\link{popgenreport}}
}
}
\details{
Check the help pages of \code{\link{popgenreport}} how to include coordinates to a genind object. The coordinates need to be projected. Latlongs are not valid, because Euclidean distances are calcuated based on these coordinates. For an example how to convert latlongs into a projected format have a look at the vignette that comes with this package. The friction needs to be a raster and needs to be in the same projection as the genind object. Also the type of genetic distance to be used needs to be specified.
}
\value{
For distance matrices are returned. Pairwise Euclidean distances between subpopulations/individuals, cost distances, path lengths and genetic distances. Also following the approach of Wassermann et al. 2010 a series of partial mantel tests are performed.
}
\references{
Cushman, S., Wasserman, T., Landguth, E. and Shirk, A. (2013). Re-Evaluating Causal Modeling with Mantel Tests in Landscape Genetics. Diversity, 5(1), 51-72. 

Landguth, E. L., Cushman, S. A., Schwartz, M. K., McKelvey, K. S., Murphy, M. and Luikart, G. (2010). Quantifying the lag time to detect barriers in landscape genetics. Molecular ecology, 4179-4191. 

Wasserman, T. N., Cushman, S. A., Schwartz, M. K. and Wallin, D. O. (2010). Spatial scaling and multi-model inference in landscape genetics: Martes americana in northern  Idaho. Landscape Ecology, 25(10), 1601-1612.  
}
\author{
Bernd Gruber (bernd.gruber@canberra.edu.au)
}


\seealso{
  \code{\link{popgenreport}}, \code{\link{wassermann}},
  \code{\link{genleastcost}}
}
\examples{
##Not: run
# require(raster)
# data(landgen)
# data(fric.raster)
# landgenreport(landgen, fric.raster, gen.distance="D", mk.mantel=T)
}