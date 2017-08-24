\name{gauss_weighting}

\Rdversion{1.1}

\alias{gauss_weighting}

\title{
  Creation of a gaussian weighting filter for each element
}

\description{
  Given an entropy object (result of the \code{\link{entropy}} or of the \code{\link{rotamer_entropy}}function), creates a vector with a gaussian weighting of each element based on the entropy value. The vector will be used to give different weights to each element in the \code{\link{centered_pca}} function.
}

\usage{
  gauss_weighting(entropy, L= 0.1, lambda= 0)
}

\arguments{
 \item{entropy}{
  An object created by the \code{\link{entropy}} function
 }
 \item{L}{
  A value indicating the gaussian curve's width. It is recommanded to be between 0.05 and 0.2. (0.1 by default)
 }
 \item{lambda}{
  A value between 0 and 1 that will be added to the gaussian weighting score to avoid zero problems
 }
}

\details{
  The object returned by the \code{\link{entropy}} function contains an entropy score for each position.
  The weighting of each position is calculated as follows :
      \deqn{weighting[i]= e^{\frac{-(entropy[i]-0.5)^2}{L}}}
      
  Then the lambda constant is applied as follow :
      \deqn{weighting= lambda + (1-lambda) * weigthing}
}

\value{
  A numeric matrix that contains a weighting score for each position in the alignment based on gaussian function.
}

\author{
  Antoine GARNIER
}

\examples{
 align <- import.msf(system.file("msa/toy_align.msf", package = "Bios2cor"))

  #Creating ENTROPY object
  entropy <- entropy(align)
  
  filter <- gauss_weighting(entropy, L= 0.1, lambda= 1)
  
}

\keyword{weighting}
\keyword{entropy}
\keyword{gauss}
