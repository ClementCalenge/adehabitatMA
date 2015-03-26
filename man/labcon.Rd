\name{labcon}
\alias{labcon}
\title{Labelling Connected Features }
\description{
  This function attributes unique labels to pixels belonging to connected
  features on a map of class \code{SpatialPixelsDataFrame}.
}
\usage{
labcon(x)
}
\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame} with one
    column}
}
\value{
  Returns a matrix of class \code{asc}, of type \code{"factor"}, with a
  number of levels equals to the number of connected components
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link[sp]{SpatialPixelsDataFrame-class}} for further
  information on the class \code{SpatialPixelsDataFrame}}

\examples{

data(lynxjura)
sa <- lynxjura$map[,1]

## build an image with separate components
sa[[1]][sa[[1]]<50] <- NA

## show the connected components
image(sa)
image(labcon(sa))


}
\keyword{spatial}

