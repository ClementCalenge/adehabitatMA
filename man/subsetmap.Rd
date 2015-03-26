\name{subsetmap}
\alias{subsetmap}
\title{Storing a Part of a Map}
\description{
  \code{subsetmap} is used to store a part of any given map of class
  \code{SpatialPixelsDataFrame} into an other object.
}
\usage{
subsetmap(x, xlim = NULL, ylim = NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame}}
  \item{xlim}{numerical vector of length 2.  The x limits of the rectangle
    including the new map}
  \item{ylim}{numerical vector of length 2.  The y limits of the rectangle
    including the new map}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  If \code{xlim} or \code{ylim} are not provided, the function asks the
  user to click on the map to delimit the lower left corner and the
  higher right corner of the new map (see Examples).
}
\value{
  Returns an object of class \code{SpatialPixelsDataFrame}
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr},
  improvements by Jon Olav Vik}
\examples{

data(lynxjura)
map <- lynxjura$map

## limits of the new map:
xl <- c(839938.7, 858990.8)
yl <- c(2149019, 2168761)

## computation of the new map:
su <- subsetmap(map, xlim = xl, ylim = yl)
su

## Display
opar <- par(mar = c(0,0,0,0))
layout(matrix(c(1,1,1,1,1,1,1,1,2), byrow = TRUE, ncol = 3))
image(map, axes = FALSE)
polygon(c(xl[1], xl[2], xl[2], xl[1]),
        c(yl[1], yl[1], yl[2], yl[2]))
image(su, axes = FALSE)
box()

par(opar)
par(mfrow = c(1,1))


\dontrun{
 ## Interactive example 
 su <- subsetmap(map)

 image(su)
}

}
\keyword{spatial}
