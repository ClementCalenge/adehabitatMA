\name{ascgen}
\alias{ascgen}
\title{Making Raster Maps From SpatialPoints Objects}
\description{
  \code{ascgen} creates an object of class \code{SpatialPixelsDataFrame}
  using a set of points contained in an object of class
  \code{SpatialPointsDataFrame}.
}
\usage{
ascgen(xy, cellsize = NULL, nrcol = NULL, count = TRUE)
}

\arguments{
  \item{xy}{an object inheriting the class \code{SpatialPoints}}
  \item{cellsize}{the cell size (resolution) of the object of class
    \code{SpatialPixelsDataFrame} to be built}
  \item{nrcol}{the size of the square raster map to
    be created (number of rows/columns)}
  \item{count}{logical.  If \code{TRUE}, the resulting object contains
    the number of points in each cell.  If \code{FALSE}, all the cells
    are set to zero} 
}
\value{
  Returns an object of class \code{SpatialPixelsDataFrame}.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link[sp]{SpatialPixelsDataFrame-class}} for additional
  information on objects of class \code{SpatialPixelsDataFrame} and
  \code{\link[sp]{SpatialPoints-class}} for additional information on
  objects of class \code{SpatialPoints}.}

\examples{

## generates a random sample of points
xy <- matrix(runif(1000), ncol=2)

## coerce them to SpatialPoints
xy <- SpatialPoints(xy)
plot(xy)

## generate a SpatialPixelsDataFrame
## (and count the number of points)
spd <- ascgen(xy, cellsize=0.1)
image(spd)


}
\keyword{spatial}
