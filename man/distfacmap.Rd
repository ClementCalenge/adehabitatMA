\name{distfacmap}
\alias{distfacmap}
\title{Compute distances to the different levels of a factor map }
\description{
  This function computes maps of distances to patches belonging to the
  different levels of a factor variable (mapped in an object of class
  \code{SpatialPixelsDataFrame}).
}
\usage{
distfacmap(x, lev = NULL)
}
\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame} with one
    column (considered as a factor by the function)}
  \item{lev}{a vector of character strings giving the labels of the
    levels of the factor.}
}
\value{
  An object of class \code{SpatialPixelsDataFrame}.
}
\author{  Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{
  \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.  
}
\examples{

## Load meuse.grid data set and converts it to
## SpatialPixelsDataFrame
data(meuse.grid)
m <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                            data = meuse.grid)

## look at the soil type
image(m[,6])

## compute the distance to each soil type
sor <- distfacmap(m[,6], lev = c("type1","type2", "type3"))

## The results
sor
mimage(sor)

}
\keyword{spatial}
