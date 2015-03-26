\name{hist.SpatialPixelsDataFrame}
\alias{hist.SpatialPixelsDataFrame}
\title{Histograms of Mapped Variables}
\description{
  \code{hist.SpatialPixelsDataFrame} performs histograms of the
  variables mapped in objects of class \code{SpatialPixelsDataFrame}.
}
\usage{
\method{hist}{SpatialPixelsDataFrame}(x, type = c("h", "l", "b"),
                                      adjust = 1, col, border,
                                      lwd = 1, \dots)
}
\arguments{
  \item{x}{a raster map of class \code{SpatialPixelsDataFrame} }
  \item{type}{what type of plot should be drawn.  Possible types are:\cr
    * \code{"h"} for histograms,\cr
    * \code{"l"} for kernel density estimates (see \code{?density}).\cr
    * \code{"b"} for both histograms and kernel density estimates (see \code{?density}).\cr
    By default, \code{type = "h"} is used.  If \code{type = "l"} is used,
    the position of the mean of each distribution is indicated by dotted
    lines}
  \item{adjust}{if \code{type = "l"}, a parameter used to control the
    bandwidth of the density estimate (see \code{?density})}
  \item{col}{color for the histogram}
  \item{border}{color for the border of the histogram}
  \item{lwd}{if \code{type = "l"}, line width for the density estimate}
  \item{\dots}{further arguments passed to or from other methods}
}
\author{Mathieu Basille \email{basille@ase-research.org}}
\seealso{ \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.}
\examples{

data(lynxjura)
hist(lynxjura$map, type = "h")

hist(lynxjura$map, type = "l")

}
\keyword{hplot}

