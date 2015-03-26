\name{buffer}
\alias{buffer}
\title{Compute Buffers Regions}
\description{
  \code{buffer} computes buffers regions from spatial objects belonging 
  to the following classes:
  
  \code{SpatialPoints} (all pixels located within a specified distance
  of a point take the value one)
  
  \code{SpatialPointsDataFrame} with one column (this column is
  considered to be a factor, and the buffer is computed for each level
  of this factor)
  
  \code{SpatialLines} to compute buffers from lines.
  
  \code{SpatialPolygons} to compute buffers from polygons.
  
}
\usage{
buffer(xy, x, dist)
}

\arguments{
  \item{xy}{an object of class \code{SpatialPoints},
    \code{SpatialPointsDataFrame} with one column (a factor considered
    to be the identity of the points), \code{SpatialLines}, or
    \code{SpatialPolygons}.}
  \item{x}{an object inheriting the class \code{SpatialPixels} with the
    same attributes (resolution, size) as those desired for the output.}
  \item{dist}{a value of distance}
}
\value{
  An object of class \code{SpatialPixelsDataFrame}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link[sp]{SpatialPixelsDataFrame-class}} for additionnal
  information on objects of class \code{SpatialPixelsDataFrame}.
}
\examples{

data(lynxjura)

# locs is the SpatialPointsDataFrame containing the
# locations of lynx indices in the Jura mountains
locs <- lynxjura$locs
head(locs)

## just for the sake of illustration: sample 100 points
set.seed(233)
locs <- locs[sample(1:nrow(locs), 100),]

# sa is the SpatialPixelsDataFrame object containing
# maps of the study area
sa <- lynxjura$map

# Buffer of 2000 m from all points
bu <- buffer(locs, sa, 2000)
image(bu)

# displays all the pixels of the study area within 2000 m
# of a point, for each type of indices (see ?lynxjura)
buani <- buffer(locs[,2], sa, 2000)
buani
par(mar=c(0,0,2,0))
opar<-par(mfrow=c(3,4))
lapply(1:11, function(i) {
   image(buani, i)
   title(main = names(slot(buani, "data"))[i])
})
par(opar)

## buffer from a polygon
sa2 <- sa[,4]
sa2[[1]][sa2[[1]]<5000] <- NA
image(sa2)

## gets the contour line
gc <- getcontour(sa2)
plot(gc, add=TRUE)

## a buffer of 2000 metres
image(buffer(gc, sa, 2000))
plot(gc, add=TRUE)



}
\keyword{spatial}
