% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lgrMMRR.r
\name{lgrMMRR}
\alias{lgrMMRR}
\title{Multiple Matrix Regression with Randomization analysis}
\usage{
lgrMMRR(gen.mat, cost.mats, eucl.mat = NULL, nperm = 999)
}
\arguments{
\item{gen.mat}{a genetic distance matrix (e.g. output from
\code{\link{genleastcost}}}

\item{cost.mats}{a list of cost distance matrices}

\item{eucl.mat}{pairwise Euclidean distance matrix. If not specificed
ignored}

\item{nperm}{the number of permutations}
}
\value{
a table with the results of the matrix regression analysis.
(regression coefficients and associated p-values from the permutation test
(using the pseudo-t of Legendre et al. 1994).  and also r.squared from and
associated p-value from the permutation test. F.test.

Finally also the F-statistic and p-value for overall F-test for lack of fit.
}
\description{
performs Multiple Matrix Regression with Randomization analysis This method
was implemented by Wang 2013 (MMRR function see references) and also by
Sarah Goslee in package ecodist. lgrMMRR is a simple wrapper to have a more
user friendly output.
}
\details{
Performs multiple regression on distance matrices following the methods
outlined in Legendre et al. 1994 and implemented by Wang 2013.
}
\examples{

\donttest{
data(landgen)
library(raster)
fric.raster <- readRDS(system.file("extdata","fric.raster.rdata", package="PopGenReport"))
glc <- genleastcost(landgen, fric.raster, "D", NN=4, path="leastcost")
lgrMMRR(glc$gen.mat, glc$cost.mats, glc$eucl.mat, nperm=999)
}
}
\references{
Legendre, P.; Lapointe, F. and Casgrain, P. 1994. Modeling brain
evolution from behavior: A permutational regression approach. Evolution 48:
1487-1499.

Lichstein, J. 2007. Multiple regression on distance matrices: A multivariate
spatial analysis tool. Plant Ecology 188: 117-131.

Wang,I 2013. Examining the full effects of landscape heterogeneity on
spatial genetic variation: a multiple matrix regression approach for
quantifying geographic and ecological isolation. Evolution: 67-12:
3403-3411.
}
\seealso{
\code{\link{popgenreport}},
\code{\link{genleastcost}}, \code{\link{landgenreport}},
\code{\link{wassermann}}
}
\author{
Bernd Gruber (bernd.gruber@canberra.edu.au) using the
implementation of Wang 2013.
}
