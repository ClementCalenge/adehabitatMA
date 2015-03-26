\name{lowres}
\alias{lowres}
\title{Reducing the Resolution of a Map}
\description{
  \code{lowres} is used to reduce the resolution of
  maps of class \code{SpatialPixelsDataFrame}.
}
\usage{
lowres(x, np = 2, which.fac=NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame}}
  \item{np}{a number giving the number of pixels to merge together
    (see below)}
  \item{which.fac}{a vector containing the indices of the columns of
    \code{x}, which should be considered as a factor}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  The function merges together squares of \code{np * np} pixels.  For
  variables of type \code{"numeric"}, the function averages the value of
  the variable.  For maps of type \code{"factor"}, the function gives
  the most frequent level in the square of \code{np * np} pixels.  When
  several levels are equally represented in the square of \code{np * np}
  pixels, the function randomly samples one of these levels.
}
\value{
  Returns an object of class \code{SpatialPixelsDataFrame}.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link[sp]{SpatialPixelsDataFrame-class}} for further
  information on objects of class \code{SpatialPixelsDataFrame}.}
\examples{
data(meuse.grid)
m <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                            data = meuse.grid)
m

m <- m[,3:6]

## The initial image
image(m,3)

## The transformed image
m2 <- lowres(m, np = 4)
image(m2, 3)

}
\keyword{spatial}
