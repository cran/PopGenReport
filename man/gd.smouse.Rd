\name{gd.smouse}
\alias{gd.smouse}
\title{Individual genetic distance calculation based on Smouse and Peakall 1999}
\description{
Calculate pairwise genetic distances between all individuals using the individual genetic distance measure of Smouse and Peakall (1999). This function should produce the same results as the individual genetic distances calculated using GenAlEx and choosing the interpolate missing data option. Note that depending on your computers capabilities, you may run into memory errors when this function is used on datasets with large numbers of individuals (>1000).
}
\usage{
gd.smouse(population, verbose = TRUE)
}
\arguments{
  \item{population}{this is the \code{\link{genind}} object that the analysis will be based on.}
  \item{verbose}{information on progress. For large data sets this could be informative as it allows the user to estimate the amount of time remaining until the function is complete. The counter shows the number of the population for which genetic distances are currently being determined.
}
}
\value{
Returns pairwise individual genetic distances for each individual within a population. 
}
\references{
Smouse PE, Peakall R. 1999. Spatial autocorrelation analysis of individual multiallele and multilocus genetic structure. Heredity 82: 561-573.}
\author{
Aaron Adamack, aaron.adamack@canberra.edu.au
}
\seealso{
\code{\link{popgenreport}}
}
\examples{
#not run:
#data(bilby)
#popgenreport(bilby, mk.gd.smouse = TRUE, mk.pdf=FALSE)
#to get a pdf output you need to have a running Latex version installed on your system.
#popgenreport(bilby, mk.gd.smouse = TRUE, mk.pdf=TRUE)
}
