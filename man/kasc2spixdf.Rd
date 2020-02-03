\name{kasc2spixdf}
\alias{kasc2spixdf}
\alias{asc2spixdf}
\alias{area2spol}
\alias{spol2area}
\alias{attpol2area}
\title{ Conversion of old classes (adehabitat) to new classes (sp,adehabitatMA)}
\description{

  These functions convert maps of classes available in adehabitat toward
  classes available in the package \code{sp} and conversely.
  
  \code{kasc2spixdf} converts an object of class \code{kasc} into an
  object of class \code{SpatialPixelsDataFrame}.
  
  \code{asc2spixdf} converts an object of class \code{asc} into an
  object of class \code{SpatialGridDataFrame}.
  
  \code{area2spol} converts an object of class \code{area} into an
  object of class \code{SpatialPolygons}.
  
  \code{spol2area} converts an object of class \code{SpatialPolygons} or
  \code{SpatialPolygonsDataFrame} into an object of class
  \code{area}.
  
  \code{attpol2area} gets the data attribute of an object of class
  \code{SpatialPolygonsDataFrame} and stores it into a data frame.
}
\usage{
kasc2spixdf(ka)
asc2spixdf(a)
area2spol(ar)
spol2area(sr)
attpol2area(srdf)
}
\arguments{
  \item{ka}{an object of class \code{kasc}.}
  \item{a}{an object of class \code{asc}.}
  \item{ar}{an object of class \code{area}.}
  \item{sr}{an object of class \code{SpatialPolygons} or
    \code{SpatialPolygonsDataFrame}.}
  \item{srdf}{an object of class \code{SpatialPolygonsDataFrame}.}
}
\details{
  We describe here more in detail the functions \code{spol2area} and
  \code{attpol2area}.  Objects of class \code{area} do not deal with
  holes in the polygons, whereas the objects of class
  \code{SpatialPolygons} do.  Therefore, when holes are present in the
  \code{SpatialPolygons} object passed as argument, the function
  \code{spol2area} ignore them and returns only the external contour of
  the polygon (though a warning is returned).
}
\author{ Clement Calenge \email{clement.calenge@ofb.gouv.fr}}
\keyword{hplot}
