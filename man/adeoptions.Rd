\name{adeoptions}
\alias{adeoptions}
\title{Setting options for the adehabitat* package }
\description{
  The function \code{adeoptions} defines the options for the package
  \code{adehabitatMA}, and more generally, for all brother packages
  (\code{adehabitatHR}, \code{adehabitatHS} and \code{adehabitatLT}).
}
\usage{
adeoptions(\dots)
}
\arguments{
  \item{\dots}{a list of named elements giving the value of options for
    the package \code{adehabitatMA}}
}
\details{
  This functions sets the value of components of the internal list
  \code{.adeoptions}, which contains the values of the options.   Two
  options are currently implemented:\cr
  
  \code{epsilon}: The smallest significant distance between two points
    above which we should consider to deal with two separate
    locations.\cr
  
  \code{shortprint}:  logical.  Defines whether special \code{show}
    methods should be used when printing objects belonging to classes of
  the package \code{sp}
}
\value{
  a list of options for the package \code{adehabitatMA}
}
\author{ Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\examples{

## load the data
data(lynxjura)

## short printing method
oldopt <- adeoptions(shortprint=TRUE)

lynxjura$map

## 'sp' print method for Spatial objects
adeoptions(shortprint=FALSE)

lynxjura$map

## original options
adeoptions(oldopt)

lynxjura$map


}
\keyword{environment}
