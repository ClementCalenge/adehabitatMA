\name{mimage}
\alias{mimage}
\title{Displaying Multi-layer Raster Maps }
\description{
  This function allows to display the whole content of an object of
  class \code{SpatialPixelsDataFrame}
}
\usage{
mimage(x, var = names(slot(x, "data")), col = gray((240:1)/256), mfrow = NULL)
}
\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame}}
  \item{var}{The names or index of the variables to be plotted }
  \item{col}{a vector of colors to be used for plotting}
  \item{mfrow}{The parameter \code{mfrow} for the resulting graph (see
    \code{help(par)} }
}
\author{  Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{par}} for information about \code{mfrow}, and
  \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.}
\examples{

data(lynxjura)

lynxjura$map

mimage(lynxjura$map)

mimage(lynxjura$map, c("forets", "routes"),
       col=grey(seq(0,1, length=100)))

}
\keyword{spatial}

