\name{gd_kosman}
\alias{gd_kosman}
\title{Individual genetic distance calculation based on Kosman & Leonhard 2005}
\description{
 Calculates pairwise genetic distances between all individuals using the individual genetic distance measure of Kosman and Leonard (2005). This function is similiar to the dist.codom in the package mmod. The two functions differ in their treatment of individuals with missing data. dist.codom omits individuals from the calculation of pairwise individual genetic distances while this function includes individuals with missing data. This is done by simply calculating the mean individual pairwise genetic distance over all loci for which there are values. Note that depending on your computers capabilities, you may run into memory errors when this function is used on datasets with large numbers of individuals (>1000).
}
\usage{
gd_kosman(population, verbose = TRUE)
}
\arguments{
  \item{population}{this is the \code{\link{genind}} object the analysis will be based on.}
  \item{verbose}{ information on progress. For large data sets this could be informative as it allows to estimate the required time until completion.}
}
\value{
Returns pairwise individual genetic distances for each individual within a population.
}
\references{
Kosman E., Leonard K.J. 2005. Similarity coefficients for molecular markers in studies of genetic relationships between individuals for haploid, diploid, and polyploidy species. Molecular Ecology 14:415-424.
}
\author{
Aaron Adamack, aaron.adamack@canberra.edu.au

}
\seealso{
\code{\link{popgenreport}}
}
\examples{
data(bilby)
#here we use only the first 50 individuals to speep up the example
popgenreport(bilby[1:50], mk.gd.kosman = TRUE, mk.pdf=FALSE)

#to get a pdf output you need to have a running Latex version installed on your system.
#popgenreport(bilby, mk.gd.kosman = TRUE, mk.pdf=TRUE)


}
