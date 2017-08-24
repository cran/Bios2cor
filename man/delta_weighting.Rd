\name{delta_weighting}

\Rdversion{1.1}

\alias{delta_weighting}

\title{
  Creation of a Delta weighting filter for each element
}

\description{
  Given an entropy object (result of the \code{\link{entropy}} or of the \code{\link{rotamer_entropy}}function), creates a vector with a delta weighting of each element based on the entropy value. The vector will be used to give different weights to each element in the \code{\link{centered_pca}} function.
}

\usage{
  delta_weighting(entropy, Smin= 0.4, Smax= 0.6)
}

\arguments{
 \item{entropy}{
  An object created by the \code{\link{entropy}} function
 }
 \item{Smin}{
  A value indicating the minimum entropy value. (Smin= 0.4 by default)
 }
 \item{Smax}{
  A value indicating the maximum entropy value. (Smax= 0.6 by default)
 }
}

\details{
  The object returned by the \code{\link{entropy}} or the \code{\link{rotamer_entropy}} function contains an entropy score for each element.
  The weighting of each element is calculated as follow :
      \deqn{ weighting[i] = \left\{ \begin{array}{ll}  1, & Smin < entropy[i] < Smax \\  0, & otherwise \end{array}  \right. }
}

\value{
  A vecto that contains a weighting score for each element (position in sequence alignment or side chain dihedral angle in trajectory).
}

\author{
  Antoine GARNIER
}

\examples{
  align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating ENTROPY object
  entropy <- entropy(align)
  
  filter <- delta_weighting(entropy, Smin= 0.4, Smax= 0.6)
  
}

\keyword{weighting}
\keyword{entropy}
\keyword{delta}
