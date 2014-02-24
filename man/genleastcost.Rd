\name{genleastcost}
\alias{genleastcost}
\title{
Least-cost path analysis based on a friction matrix
}
\description{
This function calculates the pairwise distances (Euclidean, Least cost path and genetic distances) of populations using a friction matrix and a spatial genind object.
The genind object needs to have coordinates in the same projected coordinate system
}
\usage{
genleastcost(cats, fr.raster, gen.dist)
}
\arguments{
  \item{cats}{a spatial genind object. see ?popgenreport how to provide coordinates in genind objects
}
  \item{fr.raster}{a friction matrix
}
  \item{gen.dist}{specification which genetic distance method should be used to calculate pairwise genetic distances between population ( "D", "Gst.Nei", "Gst.Hedrick", "Smouse" or "Kosman")
}
}
\details{
to be written
}
\value{
returns a list that consists of four pairwise distance matrixes (Euclidean, Cost, length of path and genetic) and the actual paths as spatial line objects
}
\references{
Cushman et al.
Landguth et al.
}
\author{
Bernd Gruber
}
\seealso{
 \code{\link{landgenreport}}, \code{\link{popgenreport}}
}
\examples{
##Not: run
# require(raster)
# data(landgen)
# data(fric.raster)
# glc <- genleastcost(landgen, fric.raster, "D")
# wassermann(eucl.mat = glc$eucl.mat, cost.mat = glc$cost.mat, gen.mat = glc$gen.mat)
}
