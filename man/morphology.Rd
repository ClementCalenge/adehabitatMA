\name{morphology}
\alias{morphology}
\title{Morphology: Erosion or Dilatation of Features on a Raster Map}
\description{
  \code{morphology} performs morphological operations on images
  of class \code{SpatialPixelsDataFrame}.
}
\usage{
morphology(x, operation = c("erode", "dilate"), nt = 5)
}
\arguments{
  \item{x}{a matrix of class \code{SpatialPixelsDataFrame} with one
    column}
  \item{operation}{a character string indicating the operation to be
    processed: either \code{"erode"} or \code{"dilate"}}
  \item{nt}{the number of times that the operation should be processed}
}
\value{
  Returns an object of class \code{SpatialPixelsDataFrame} with one
  column, containing 1 when the pixel belong to one feature of the image
  and \code{NA} otherwise (see examples).
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link[sp]{SpatialPixelsDataFrame-class}} for further
  information on objects of class \code{SpatialPixelsDataFrame}.} 
\examples{

data(meuse.grid)
a <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                            data = meuse.grid)


## dilatation
toto1 <- morphology(a, operation = "dilate", nt = 1)
toto2 <- morphology(a, operation = "dilate", nt = 2)
toto3 <- morphology(a, operation = "dilate", nt = 3)
toto5 <- morphology(a, operation = "dilate", nt = 5)
colo <- grey((1:5)/6)
image(toto5, col = colo[1])
image(toto3, add = TRUE, col = colo[2])
image(toto2, add = TRUE, col = colo[3])
image(toto1, add = TRUE, col = colo[4])
image(a, add = TRUE)

## erosion
colo <- grey((1:12)/13)
image(a, col = "black")
for (i in 1:12) {
  toto <- morphology(a, operation = "erode", nt = i)
  image(toto, add = TRUE, col = colo[i])
}


}
\keyword{spatial}

