\name{count.points}
\alias{count.points}
\title{Number of Points in Each Pixel of a Raster Map}
\description{
  \code{count.points} counts the number of points in
  each pixel of a raster map inheriting the class
  \code{SpatialPixels}.
}
\usage{
count.points(xy, w)
}
\arguments{
  \item{xy}{an object of class \code{SpatialPoints}, or
    \code{SpatialPointsDataFrame} with one column.  In the latter case,
    the column is considered as a factor giving, for each point, the
    membership of the point to a set.}
  \item{w}{an object inheriting the class \code{SpatialPixels}.}
}
\value{
  an object of class \code{SpatialPixelsDataFrame} containing the number
  of points in each cell of the raster map.  If \code{xy} is a
  \code{SpatialPointsDataFrame} with one column (a factor), the
  resulting object contains one column per level of this factor.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.}

\examples{

data(lynxjura)

# locs is the SpatialPointsDataFrame containing the
# locations of presence indices of Lynx in the Jura mountains (France)
locs <- lynxjura$locs
head(locs)

# sa is the SpatialPixelsDataFrame object containing
# maps of the study area
sa <- lynxjura$map

# Count all points
cp <- count.points(locs, sa)
cp
image(cp)

# count the points per type of lynx presence indices:
cp <- count.points(locs[,2], sa)
cp


}
\keyword{spatial}
