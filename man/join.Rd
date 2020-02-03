\name{join}
\alias{join}
\title{Finds the Value of Mapped Variables at some Specified Locations
  (Spatial Join)}
\description{
  \code{join} finds the value of a mapped variable at some specified
  locations.
}
\usage{
 join(xy, x) 
}

\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame}}
  \item{xy}{an object of class \code{SpatialPointsDataFrame}}
}
\value{
  If only one variable is mapped in \code{x}, a vector with length
  equals to the number of points in \code{xy}.

  If only several variables are mapped in \code{x}, a \code{data.frame}
  with a number of columns equal to the number of variables in the
  object of class \code{SpatialPixelsDataFrame}, and with each row
  corresponding to the rows of \code{xy}. 
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr} }

\seealso{ \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.
  \code{\link[sp]{SpatialPoints-class}} for additionnal information on
  objects of class \code{SpatialPoints}.} 

\examples{
data(lynxjura)

## show the data
image(lynxjura$map)
points(lynxjura$locs)

## join the data to the maps:
res <- join(lynxjura$locs, lynxjura$map)
head(res)

}
\keyword{programming}
\keyword{hplot}
