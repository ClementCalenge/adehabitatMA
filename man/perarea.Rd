\name{calcperimeter}
\alias{calcperimeter}
\title{Compute Perimeters of Objects of Class "SpatialPolygonsDataFrame"
  and "PolyLinesDataFrame"}
\description{
  \code{calcperimeter} computes the perimeters of polygons in objects of
  class \code{SpatialPolygonsDataFrame} and \code{PolyLinesDataFrame}.
}
\usage{
  calcperimeter(x)
}
\arguments{
  \item{x}{object of class \code{SpatialPolygonsDataFrame} or
    \code{PolyLinesDataFrame}}
}
\value{
  an object of class \code{SpatialPolygonsDataFrame} or
  \code{PolyLinesDataFrame} with an additionnal column containing the
  perimeter of the polygons/polylines.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link[sp]{SpatialPolygonsDataFrame-class}} for
  additionnal information on objects of class
  \code{SpatialPolygonsDataFrame}.}
\examples{
data(meuse.grid)
a <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                            data = meuse.grid)


## the contour of the map
gc <- getcontour(a[,1])
plot(gc)

## Transforms the SpatialPolygons into SpatialPolygonsDataFrame
gc <- SpatialPolygonsDataFrame(gc, data.frame(x=1))

## The perimeter of this map (in units of the data):
ii <- calcperimeter(gc)
as.data.frame(ii)



}
\keyword{hplot}

