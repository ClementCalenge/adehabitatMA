\name{explore}
\alias{explore}
\title{Interactive Exploration of Maps of Class 'SpatialPixelsDataFrame'
  (requires the package tkrplot)}
\description{
  This interface allows to explore distances, values, etc. on a map of
  class \code{SpatialPixelsDataFrame}.
}
\usage{
explore(ka, coltxt="black",
        hscale = 1, vscale = 1,
        panel.last = NULL, \dots)
}
\arguments{
  \item{ka}{An object of class \code{kasc} }
  \item{coltxt}{character.  the color of the text to be printed}
  \item{hscale}{passed to tkrplot}
  \item{vscale}{passed to tkrplot}
  \item{panel.last}{an expression to be evaluated after plotting has taken
    place}
  \item{\dots}{additional parameters to be passed to the function
    \code{image}} 
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{  \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.}
\examples{
\dontrun{

data(meuse.grid)
m <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                            data = meuse.grid)
explore(m[,3:6])

data(lynxjura)
explore(lynxjura$map, panel.last=function() points(lynxjura$locs))


}
}
\keyword{hplot}
\keyword{spatial}
