\name{puechabonsp}
\docType{data}
\alias{puechabonsp}
\title{ Radio-Tracking Data of Wild Boar }
\description{
This data set stores the results of the monitoring of 4 wild boars at
Puechabon (Mediterranean habitat, South of France).
These data have been collected by Daniel Maillard (Office national de la
chasse et de la faune sauvage).
}
\usage{
data(puechabonsp)
}
\details{
  The list \code{puechabonsp} has two components:
  
  \code{puechabonsp$map}
  is an  object of class \code{SpatialPixelsDataFrame} that describes several
  variables on the study area.

  \code{puechabonsp$relocs} is an object of class \code{SpatialPointsDataFrame}
  containing the relocations of the wild boar resting sites in
  summer. Information on wild boars is provided by factors \code{Name},
  \code{Sex}, \code{Age}.
}
\references{
Maillard, D. (1996). \emph{Occupation et utilisation de la garrigue et
  du vignoble mediterraneens par le Sanglier}.
Universite d'Aix-Marseille III: PhD thesis.
}

\keyword{datasets}
