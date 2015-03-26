\name{getcontour}
\alias{getcontour}
\title{Computes the Contour Polygon of a Raster Object}
\description{
  \code{getcontour} computes the contour polygon of a raster object
  of class \code{SpatialPixelsDataFrame}.  When the object is made of
  several parts, the function returns one polygon per part.
}
\usage{
getcontour(sp)
}

\arguments{
  \item{sp}{an object of class \code{SpatialPixelsDataFrame}}
}
\value{
  Returns an object of class \code{SpatialPolygons}.
}
\section{Warning }{
  Holes in the polygons are not taken into account by the function.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.}

\examples{

data(meuse.grid)
a <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                            data = meuse.grid)


## the contour of the map
gc <- getcontour(a[,1])
plot(gc)



}
\keyword{programming}
\keyword{hplot}
