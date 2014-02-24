\name{wassermann}
\alias{wassermann}
\title{
Partial Mantel tests on costdistance matrices
}
\description{
This function implements the Causal modelling approach as suggested by Wassermann et al. 2010 and Cushman et al. 2010. It tests for the effect of landscape features using a cost distance matrix on the genetic structure of subpopulation/individuals. For a complete 
}
\usage{
wassermann(gen.mat, cost.mat, eucl.mat, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen.mat}{ pairwise genetic distance matrix
}
  \item{cost.mat}{ pairwise cost distance matrix
}
  \item{eucl.mat}{ pairwise cost distance matrix
}
  \item{plot}{ switch for control plots of the partial mantel test
}
}
\details{see \code{\link{landgenreport}}
}
\value{
A table with the results of the partial mantel test.
}
\references{
Wassermann, T.N., Cushman, S. A., Schwartz, M. K. and Wallin, D. O. (2010). Spatial scaling and multi-model inference in landscape genetics: Martes americana in northern Idaho. Landscape Ecology, 25(10), 1601-1612.
}
\author{
Bernd Gruber (bernd.gruber@canberra.edu.au)
}

\seealso{
  \code{\link{popgenreport}}, \code{\link{genleastcost}}, \code{\link{landgenreport}}
}
\examples{
##Not: run
# require(raster)
# data(landgen)
# data(fric.raster)
# glc <- genleastcost(landgen, fric.raster, "D")
# wassermann(eucl.mat = glc$eucl.mat, cost.mat = glc$cost.mat, gen.mat = glc$gen.mat)
}
